//! This implements formula 3+4 of the document.
use std::mem;
use std::path::Path;
use std::thread;

use anyhow::Result;
use bio::stats::LogProb;

use rayon::prelude::*;
use statrs::function::beta::ln_beta;
// use datetime::Instant;
use std::time::{Duration, Instant};

use rayon::ThreadPoolBuilder;
use std::fs;
use crate::errors::Error;
use crate::preprocess::LnBetaCache;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::compute_grid;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;

pub(crate) fn sample_expression(
    preprocessing: &Path,
    sample_id: &str,
    epsilon: LogProb,
    c: f64,
    threads: usize,
    out_dir_path: &Path,
) -> Result<()> {
    let time0 = Instant::now();
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let sample_ids = preprocessing
        .scale_factors()
        .keys()
        .cloned()
        .collect::<Vec<_>>();

    let mean_disp_estimates =
        preprocessing
            .mean_disp_estimates()
            .get(sample_id)
            .ok_or(Error::UnknownSampleId {
                sample_id: sample_id.to_owned(),
            })?;

    let s_j = *preprocessing.scale_factors().get(sample_id).unwrap();
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let base_db_path = out_dir_path.to_str().unwrap().to_string();

    let query_points_per_feature = query_points::QueryPointsPerFeature::new(&preprocessing, c);

    // -------------------------------
    // 1. HPC SAFEGUARD
    // -------------------------------
    let max_threads = std::cmp::min(threads, 32);
    let custom_pool = ThreadPoolBuilder::new()
        .num_threads(max_threads)
        .build()?;


    // -------------------------------
    // 2. SCATTER PHASE
    // -------------------------------
    let temp_paths: Vec<String> = (0..max_threads)
        .map(|i| format!("{}_temp_{}.duckdb", base_db_path, i))
        .collect();
    print!("Time taken for initialization: {:?}\n", time0.elapsed());


    custom_pool.install(|| {
        feature_ids
            // .into_par_iter()
            .par_chunks(10)
            .map_init(
                 // INIT → runs once per thread
            {
                // let time1 = Instant::now();
                let temp_paths = temp_paths.clone();

                move || {
                    let thread_idx = rayon::current_thread_index().unwrap();

                    let path = &temp_paths[thread_idx];

                    let conn = ProbDistribution2d::open(path)
                        .expect("Failed to open DuckDB");

                    (path.clone(), conn)
                }

            },

            // WORK
            // |(temp_path, conn), (i, feature_id)| {
            |(temp_path, conn), chunk| {
                // collect all results for this chunk
                let mut time_chunk = Instant::now();
                let mut batch_results: Vec<(String, Vec<(f64, f64, f64)>)> = Vec::with_capacity(chunk.len());

                for &(i, ref feature_id) in chunk {
                    let mut time1 = Instant::now();
                    let d_ij = mean_disp_estimates.means()[i];

                    let t_ij = if let Some(t) = mean_disp_estimates.dispersions()[i] {
                        t
                    } else if let Some(t) = preprocessing.interpolate_dispersion(i) {
                        t
                    } else {
                        println!("skipped {:?}", feature_id);
                        continue;
                    };

                        let query_points = query_points_per_feature.get(i);
                        let mu_ik_points = query_points.all_mu_ik();
                        let start_points_theta_i = query_points.thetas();

                        // Pre-calculate loop invariants for the "left" nb distribution
                        let scaled_d_ij = d_ij / s_j;
                        let n_left = 1.0 / t_ij;
                        let b_left = ln_beta(scaled_d_ij + 1.0, n_left);
                        let left_const_term = (scaled_d_ij + n_left).ln();

                        let calc_prob = |m, theta_i, theta_idx| {
                            likelihood_mu_ik_theta_i(
                                scaled_d_ij,
                                m,
                                t_ij,
                                theta_i,
                                theta_idx,
                                epsilon,
                                n_left,
                                b_left,
                                left_const_term,
                                &preprocessing,
                            )
                        };
                        let probs =
                            compute_grid(&mu_ik_points, &start_points_theta_i, calc_prob);
                        print!("Calculation time taken for feature {}: {:?}\n", feature_id, time1.elapsed());
                        time1 = Instant::now();

                        batch_results.push((feature_id.to_string(), probs));
                }

                let write_start = Instant::now();

                // ONE writer per chunk
                let mut writer = ProbDistribution2d::with_connection(
                    &*conn,
                    "batch", // optional grouping key
                )
                .expect("writer");

                writer.write_batch(&batch_results).expect("write failed");
                println!(
                    "Chunk write time: {:?} (total chunk {:?})\n",
                    write_start.elapsed(),
                    time_chunk.elapsed()
                );
            },
        )
        .count(); // Force execution of the parallel iterator
    });

    let time2 = Instant::now();
    print!("Time taken for scatter phase: {:?}\n", time2 - time0);
    // -------------------------------
    // 3. GATHER PHASE
    // -------------------------------
    let final_conn = duckdb::Connection::open(&base_db_path)?;
    ProbDistribution2d::init_schema(&final_conn)?;

    for temp_path in &temp_paths {
        final_conn.execute(
            &format!("ATTACH '{}' AS temp_db", temp_path),
            [],
        )?;

        final_conn.execute(
            "INSERT INTO distributions SELECT * FROM temp_db.distributions",
            [],
        )?;

        final_conn.execute("DETACH temp_db", [])?;
    }

    // cleanup phase (parallel)
    temp_paths.par_iter().for_each(|path| {
        let _ = fs::remove_file(path);
    });
    print!("Time taken for gather phase: {:?}\n", Instant::now() - time2);

    print!("Total time taken: {:?}\n", time0.elapsed());
    Ok(())
}


/// Inner of equation 3/4 in the document.
fn likelihood_mu_ik_theta_i(
    scaled_d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    theta_idx: usize,
    _: LogProb, // epsilon
    n_left: f64,
    b_left: f64,
    left_const_term: f64,
    preprocessing: &Preprocessing,
) -> LogProb {
    if mu_ik == 0. {
        return LogProb::ln_zero();
    }
    let mut max_prob = LogProb::ln_zero();
    let mut probs = Vec::with_capacity(2048);
    let cache: &LnBetaCache = &preprocessing.ln_beta_caches()[theta_idx];
    let threshold = LogProb(0.001_f64.ln());

    let nb_right = NegBinomPrepared::new(mu_ik, theta_i, cache);

    // fast inline closure for the left pmf
    let calc_left_pmf = |x_f64: f64| -> LogProb {
        let p = n_left / (n_left + x_f64);
        let mut p1 = if n_left > 0.0 { n_left * p.ln() } else { 0.0 };
        let mut p2 = if scaled_d_ij > 0.0 { scaled_d_ij * (1.0 - p).ln() } else { 0.0 };

        if p1 < p2 {
            std::mem::swap(&mut p1, &mut p2);
        }
        LogProb((p1 - b_left + p2) - left_const_term)
    };


    for x in 0..200 {
        // let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        // let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
        let calced_prob = calc_left_pmf(x as f64) + nb_right.ln_pmf(x as f64);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        probs.push(calced_prob);
    }
    for x in 200..10000 {
        // let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        // let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
        let calced_prob = calc_left_pmf(x as f64) + nb_right.ln_pmf(x as f64);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        if calced_prob - max_prob < threshold {
            let result = LogProb::ln_sum_exp(&probs);
            return result;
        }
        probs.push(calced_prob);
    }
    let nb_right = NegBinomPreparedUncached::new(mu_ik, theta_i);
    for x in 10000.. {
        // let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        // let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
        let calced_prob = calc_left_pmf(x as f64) + nb_right.ln_pmf(x as f64);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        if calced_prob - max_prob < threshold {
            break;
        }
        probs.push(calced_prob);
    }
    let result = LogProb::ln_sum_exp(&probs);
    return result;
}

pub(crate) fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);

    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    LogProb((p1 - b + p2) - (x + n).ln())
}

struct NegBinomPreparedUncached {
    n: f64,
    ln_1mp: f64,
    p1: f64,
}

impl NegBinomPreparedUncached {
    fn new(mu: f64, theta: f64) -> Self {
        let n = 1.0 / theta;
        let p = n / (n + mu);
        Self {
            n,
            ln_1mp: (1.0 - p).ln(),
            p1: n * p.ln(),
        }
    }

    #[inline]
    fn ln_pmf(&self, x: f64) -> LogProb {
        let b = ln_beta(x + 1.0, self.n);
        let mut p1 = self.p1;
        let mut p2 = x * self.ln_1mp;

        if p1 < p2 {
            mem::swap(&mut p1, &mut p2);
        }
        LogProb((p1 - b + p2) - (x + self.n).ln())
    }
}

struct NegBinomPrepared<'a> {
    n: f64,
    ln_1mp: f64,
    p1: f64,
    beta_cache: &'a LnBetaCache,
}

impl<'a> NegBinomPrepared<'a> {
    fn new(mu: f64, theta: f64, beta_cache: &'a LnBetaCache) -> Self {
        let n = 1.0 / theta;
        let p = n / (n + mu);
        Self {
            n,
            ln_1mp: (1.0 - p).ln(),
            p1: n * p.ln(),
            beta_cache,
        }
    }

    #[inline]
    fn ln_pmf(&self, x: f64) -> LogProb {
        // let b = ln_beta(x + 1.0, self.n);
        let b = self.beta_cache.get(x as usize);
        let mut p1 = self.p1;
        let mut p2 = x * self.ln_1mp;

        if p1 < p2 {
            mem::swap(&mut p1, &mut p2);
        }
        LogProb((p1 - b + p2) - (x + self.n).ln())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn test_neg_binom() {
        assert_relative_eq!(
            neg_binom(0., 10., 2.38).exp(),
            0.2594752460369642,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            neg_binom(0., 30., 2.38).exp(),
            0.16542351363026533,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            neg_binom(0., 30., 500.38).exp(),
            0.9809648435381609,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_neg_binom_prepared_matches_old() {
        let test_cases = [
            (0.0, 10.0, 2.38),
            (0.0, 30.0, 2.38),
            (0.0, 30.0, 500.38),
            (5.0, 10.0, 0.5),
            (20.0, 50.0, 1.2),
            (200.0, 50.0, 50.0),
            (200.0, 150.0, 1.2),
            (200.0, 350.0, 1.2),
        ];

        for (x, mu, theta) in test_cases {
            let old = neg_binom(x, mu, theta).exp();

            let nb = NegBinomPrepared::new(mu, theta);
            let new = nb.ln_pmf(x).exp();

            assert_relative_eq!(old, new, epsilon = 1e-14);
        }
    }
}

// def neg_binom(x, mu, theta):
//     """ Own implementation of the negative negative binomial distribution using betaln
//     """
//     n = 1. / theta
//     p = n / (n + mu)
//     p1 = n*np.log(p) if n > 0 else 0
//     p2 = x*np.log(1-p) if x > 0 else 0
//     b = betaln(x + 1, n)
//     if (p1 < p2):
//         return exp(p2 - b + p1) / (x+n)
//     else:
//         return exp(p1 - b + p2) / (x+n)
