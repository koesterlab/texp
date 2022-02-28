//! This implements formula 3+4 of the document.

use std::fs;
use std::mem;
use std::path::Path;

use anyhow::Result;
use bio::stats::LogProb;
use getset::Getters;
// use itertools_num::linspace;
use noisy_float::types::N32;
use rayon::prelude::*;
use rmp_serde::Deserializer;
use serde::Deserialize as SerdeDeserialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;

use crate::common::{window, window_x, ProbDistribution, ProbDistribution2d};
use crate::common::{MeanDispersionPair, Outdir};
use crate::errors::Error;
use crate::preprocess::Preprocessing;

pub(crate) fn sample_expression(
    preprocessing: &Path,
    sample_id: &str,
    epsilon: LogProb,
    out_dir: &Path,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let mean_disp_estimates =
        preprocessing
            .mean_disp_estimates()
            .get(sample_id)
            .ok_or(Error::UnknownSampleId {
                sample_id: sample_id.to_owned(),
            })?;
    let group_means = 
        preprocessing
            .group_means();
            
    let s_j = preprocessing
        .scale_factors()
        .get(sample_id)
        .unwrap()
        .clone();

    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let out_dir = Outdir::create(out_dir)?;

    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            let debug_print = false; //**feature_id == String::from("ENST00000643797.1");
            let d_ij = mean_disp_estimates.means()[*i]; //TODO Do we need group mean mu_ik instead of sample mean
            let d_ik = group_means[*i];
            println!("i {:?}, d_ij {:?}, d_ik{:?}", i, d_ij, d_ik);
            // METHOD: If the per-sample dispersion is unknown, fall back to a mean interpolated from the other samples.
            let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                // TODO log message
                return Ok(());
            };

            let max_prob = prob_mu_ik_theta_i_x(d_ij, d_ij, d_ij, t_ij, prior.mean(), s_j);
            let prob_threshold = LogProb(*max_prob - 10.0f64.ln());

            let (mu_ik_left_window, mu_ik_right_window) = window(d_ik);


            // let mut points = LinkedList::<[f64;2]>::new();
            let mut likelihoods = ProbDistribution2d::new();

            let mut likelihood_mu_ik = |mu_ik| {
                let mut max_prob = LogProb::ln_zero();
                let mut process_window = |window: &Vec<f64>| {
                    for theta_i in window {
                        let prob =
                            likelihood_mu_ik_theta_i(d_ij, mu_ik, t_ij, *theta_i, s_j, epsilon);

                        likelihoods.insert(mu_ik, *theta_i, prob);

                        if prob > max_prob {
                            max_prob = prob;
                        }

                        // if prob < prob_threshold {
                        //     // break;
                        // }
                    }
                };
                process_window(prior.left_window());
                process_window(prior.right_window());

                max_prob
            };

            for mu_ik in mu_ik_left_window {
                if likelihood_mu_ik(mu_ik) < prob_threshold {
                    // break;
                }
            }

            for mu_ik in mu_ik_right_window {
                if likelihood_mu_ik(mu_ik) < prob_threshold {
                    // break;
                }
            }

            out_dir.serialize_value(feature_id, likelihoods)?;

            Ok(())
        })?;

    out_dir.serialize_value(
        "info",
        SampleInfo {
            sample_id: sample_id.to_owned(),
        },
    )?;

    Ok(())
}

#[derive(Debug, Deserialize, Serialize, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct SampleInfo {
    sample_id: String,
}

impl SampleInfo {
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(SampleInfo::deserialize(&mut Deserializer::new(
            fs::File::open(path)?,
        ))?)
    }
}

fn prob_mu_ik_theta_i_x(
    x: f64,
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    s_j: f64,
) -> LogProb {
    LogProb((neg_binom(d_ij, x, t_ij) * neg_binom(x, mu_ik * s_j, theta_i)).ln())
}

/// Inner of equation 3/4 in the document.
fn likelihood_mu_ik_theta_i(
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    s_j: f64,
    epsilon: LogProb,
) -> LogProb {
    // TODO determine whether this is the best window given that we also have access to mu_ik here.
    let (x_left_window, x_right_window) = window_x(d_ij);

    let prob = |x| prob_mu_ik_theta_i_x(x, d_ij, mu_ik, t_ij, theta_i, s_j);

    let max_prob = prob(d_ij);
    let threshold = LogProb(*max_prob - 10.0f64.ln());
    let is_informative = |prob: &LogProb| *prob >= threshold; // TODO maybe better relative to the maximum?

    LogProb::ln_sum_exp(
        &x_left_window
            .map(&prob)
            // .take_while(&is_informative)
            .chain(x_right_window.map(&prob))//.take_while(&is_informative))
            .collect::<Vec<_>>(),
    )
}

fn neg_binom(x: f64, mu: f64, theta: f64) -> f64 {
    let n = 1.0 / theta;
    let p = n / (n + mu);
    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    (p1 - b + p2).exp() / (x + n)
    // (p1 - b + p2) - (x + n).ln() // TODO is this the correct form for returning LogProb?
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
