//! This implements formula 3+4 of the document.
use std::mem;
use std::path::Path;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

use anyhow::Result;
use bio::stats::LogProb;

// use getset::Getters;
use rayon::prelude::*;
// use rmp_serde::Deserializer;
// use serde::Deserialize as SerdeDeserialize;
// use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;

use crate::errors::Error;
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
        mpsc::SyncSender<(String, Vec<(f64, f64, f64)>)>,
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
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            let d_ij = mean_disp_estimates.means()[*i];
            let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                println!("skipped {:?}", feature_id);
                return Ok(());
            };

            let query_points = query_points_per_feature.get(*i);

            let calc_prob = |m, t| likelihood_mu_ik_theta_i(d_ij, m, t_ij, t, s_j, epsilon);
            let mu_ik_points = query_points.all_mu_ik();
            let start_points_theta_i = query_points.thetas();

            // Compute grid in memory
            let probs = compute_grid(&mu_ik_points, &start_points_theta_i, calc_prob);

            // Send results to writer
            tx.send((feature_id.to_string(), probs)).unwrap();

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
    s_j: f64,
    _: LogProb, // epsilon
) -> LogProb {
    if d_ij != 0. && mu_ik == 0. {
        return LogProb::ln_zero();
    }
    if mu_ik == 0. {
        return LogProb::ln_zero();
    }
    let mut max_prob = LogProb::ln_zero();
    let mut probs = Vec::with_capacity(300);
    let mut x: f64 = 0.;

    let nb_right = NegBinomPrepared::new(mu_ik, theta_i);

    loop {
        let nb_left = NegBinomPrepared::new(x, t_ij);
        let calced_prob = nb_left.ln_pmf(d_ij / s_j) + nb_right.ln_pmf(x);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        probs.push(calced_prob);

        if x > 200. && calced_prob - max_prob < LogProb(0.001_f64.ln()) {
            break;
        }
        x = x + 1.;
    }
    let result = LogProb::ln_sum_exp(&probs);
    result
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

struct NegBinomPrepared {
    n: f64,
    ln_p: f64,
    ln_1mp: f64,
}

impl NegBinomPrepared {
    fn new(mu: f64, theta: f64) -> Self {
        let n = 1.0 / theta;
        let p = n / (n + mu);
        Self {
            n,
            ln_p: p.ln(),
            ln_1mp: (1.0 - p).ln(),
        }
    }

    #[inline]
    fn ln_pmf(&self, x: f64) -> LogProb {
        let b = ln_beta(x + 1.0, self.n);
        let mut p1 = self.n * self.ln_p;
        let mut p2 = if x > 0.0 { x * self.ln_1mp } else { 0.0 };

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
