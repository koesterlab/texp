//! This implements formula 3+4 of the document.
use std::mem;
use std::path::Path;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread;

use anyhow::Result;
use bio::stats::LogProb;

// use getset::Getters;
use rayon::prelude::*;
// use rmp_serde::Deserializer;
// use serde::Deserialize as SerdeDeserialize;
// use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;

use rayon::ThreadPoolBuilder;
use std::sync::atomic::{AtomicUsize, Ordering};
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
    out_dir_path: &Path,
) -> Result<()> {
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

    let db_path = out_dir_path.to_str().unwrap().to_string();

    // --- Channel setup ---
    let (tx, rx): (
        SyncSender<(String, Vec<(f64, f64, f64)>)>,
        Receiver<(String, Vec<(f64, f64, f64)>)>,
    ) = mpsc::sync_channel(10); //TODO Determine optimal buffer size (tradeoff between memory usage and writer blocking)

    // --- Spawn the writer thread ---
    let writer_handle = thread::spawn(move || {
        let conn = ProbDistribution2d::open(&db_path).unwrap();
        while let Ok((feature_id, grid)) = rx.recv() {
            let mut writer = ProbDistribution2d::with_connection(&conn, &feature_id).unwrap();
            writer.write_output(&grid);
        }
    });

    let query_points_per_feature = query_points::QueryPointsPerFeature::new(&preprocessing, c);

    // --- Parallel workers ---
    feature_ids
        // .par_iter()
        // .try_for_each(|(i, feature_id)| -> Result<()> {
        .par_chunks(10) // Process features in chunks to reduce overhead of thread spawning and channel communication
        .try_for_each(|chunk| -> Result<()> {
            for (i, feature_id) in chunk {
                let d_ij = mean_disp_estimates.means()[*i];
                let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                    t_ij
                } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                    t_ij
                } else {
                    println!("skipped {:?}", feature_id);
                    return Ok(());
                };

                let preprocessing = &preprocessing;
                let query_points = query_points_per_feature.get(*i);
                let mu_ik_points = query_points.all_mu_ik();
                let start_points_theta_i = query_points.thetas();

                let calc_prob = |m, theta_i, theta_idx| {
                    likelihood_mu_ik_theta_i(
                        d_ij,
                        m,
                        t_ij,
                        theta_i,
                        theta_idx,
                        s_j,
                        epsilon,
                        preprocessing,
                    )
                };

                // Compute grid in memory
                let probs = compute_grid(&mu_ik_points, &start_points_theta_i, calc_prob);

                // Send results to writer
                tx.send((feature_id.to_string(), probs)).unwrap();
            }


            Ok(())

        })?;

    drop(tx); // close channel
    writer_handle.join().unwrap();

    Ok(())
}

// #[derive(Debug, Deserialize, Serialize, Getters)]
// #[getset(get = "pub(crate)")]
// pub(crate) struct SampleInfo {
//     sample_id: String,
// }

// impl SampleInfo {
//     #[allow(unused)]
//     pub(crate) fn from_path(path: &Path) -> Result<Self> {
//         Ok(SampleInfo::deserialize(&mut Deserializer::new(
//             fs::File::open(path)?,
//         ))?)
//     }
// }

/// Inner of equation 3/4 in the document.
fn likelihood_mu_ik_theta_i(
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    theta_idx: usize,
    s_j: f64,
    _: LogProb, // epsilon
    preprocessing: &Preprocessing,
) -> LogProb {
    if d_ij != 0. && mu_ik == 0. {
        return LogProb::ln_zero();
    }
    if mu_ik == 0. {
        return LogProb::ln_zero();
    }
    let mut max_prob = LogProb::ln_zero();
    let mut probs = Vec::with_capacity(300);
    let cache: &LnBetaCache = &preprocessing.ln_beta_caches()[theta_idx];
    let threshold = LogProb(0.001_f64.ln());

    let nb_right = NegBinomPrepared::new(mu_ik, theta_i, cache);
    let scaled_d_ij = d_ij / s_j;

    for x in 0..200 {
        let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        probs.push(calced_prob);
    }
    for x in 200..10000 {
        let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
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
        let nb_left = NegBinomPreparedUncached::new(x as f64, t_ij);
        let calced_prob = nb_left.ln_pmf(scaled_d_ij) + nb_right.ln_pmf(x as f64);
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
