//! This implements formula 3+4 of the document.
use std::fs;
use std::mem;
use std::path::Path;

use anyhow::Result;
use bio::stats::LogProb;

use getset::Getters;
use rayon::prelude::*;
use rmp_serde::Deserializer;
use serde::Deserialize as SerdeDeserialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;
use duckdb::{Connection, params};
use std::sync::{Arc, Mutex};

use crate::common::Outdir;
use crate::errors::Error;
use crate::preprocess::Preprocessing;
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
    // let prior = preprocessing.prior()?;
    let mean_disp_estimates =
        preprocessing
            .mean_disp_estimates()
            .get(sample_id)
            .ok_or(Error::UnknownSampleId {
                sample_id: sample_id.to_owned(),
            })?;

    let s_j = preprocessing
        .scale_factors()
        .get(sample_id)
        .unwrap()
        .clone();
    println!("s_j {:?}", s_j);

    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();
    // println!("{:?} features", feature_ids.len());

    println!("out_dir_path {:?}", out_dir_path);
    let db_path = out_dir_path.to_str().unwrap();//format!("{}.duckdb", out_dir_path.to_str().unwrap());
    println!("db_path {:?}", db_path);
    // Open DuckDB file once
    let conn = Connection::open(db_path)?;
    ProbDistribution2d::init_schema(&conn)?; // ensure schema exists
    // Wrap in Arc<Mutex<Connection>> for parallel use
    let conn = Arc::new(Mutex::new(conn));

    feature_ids
        .par_iter().take(5)
        .try_for_each(|(i, feature_id)| -> Result<()> {
            // print start time of feature
            let time1 = std::time::SystemTime::now();
            println!("feature {:?} {:?} started at {:?}", i, feature_id, time1 );

            let query_points = query_points::calc_query_points(
                c,
                preprocessing.mean_disp_estimates().clone(),
                sample_ids.clone(),
                preprocessing.feature_ids().clone(),
                *i,
            );

            let d_ij = mean_disp_estimates.means()[*i]; //TODO Do we need group mean mu_ik instead of sample mean
                                                        // println!("{:?} d_ij {:?}", feature_id, d_ij / s_j);
                                                        // METHOD: If the per-sample dispersion is unknown, fall back to a mean interpolated from the other samples.
            let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                println!("skipped {:?}", feature_id);
                // TODO log message
                return Ok(());
            };
            println!("d_ij {:?}, t_ij {:?}", d_ij, t_ij);

            let calc_prob = |m, t| {
                let prob = likelihood_mu_ik_theta_i(
                    d_ij, m, // mu_ik
                    t_ij, t, // theta_i
                    s_j, epsilon,
                );
                prob
            };
            let mu_ik_points = query_points.all_mu_ik();
            let start_points_theta_i = query_points.thetas();

            println!("before insert_grid mu_ik_points.len() {:?}, start_points_theta_i.len() {:?}", mu_ik_points.len(), start_points_theta_i.len());
            // Each feature gets its own ProbDistribution2d handle
            let likelihoods = ProbDistribution2d::with_connection(conn.clone(), feature_id).unwrap();

            // Compute grid in memory
            let probs = likelihoods.compute_grid(&mu_ik_points, &start_points_theta_i, calc_prob);

            // Write results to DuckDB (mutex ensures serialized access)
            likelihoods.write_output(&probs).unwrap();

            let time2 = std::time::SystemTime::now();
            println!(
                "sample {} feature {:?} {:?} finished in duration {:?}",
                sample_id,
                i,
                feature_id,
                time2.duration_since(time1).unwrap()
            );

            Ok(())
        })?;

    Ok(())
}

#[derive(Debug, Deserialize, Serialize, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct SampleInfo {
    sample_id: String,
}

impl SampleInfo {
    #[allow(unused)]
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
    // println!("x {:?}, d_ij {:?}, mu_ik {:?}, t_ij {:?}, theta_i {:?}, s_j {:?}", x, d_ij, mu_ik, t_ij, theta_i, s_j);
    //METHOD /s_j?? mu_ik*s_j würde einluss von t_ij ändern, weil mu geändert wird.
    let left = neg_binom(d_ij / s_j, x, t_ij);
    let right = neg_binom(x, mu_ik, theta_i);
    // let right =  LogProb::from(Poisson::new(mu_ik).unwrap().ln_pmf(x as u64));
    let result = left + right;
    // println!("result {:?}", result);
    result
}

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
    let mut result = LogProb::ln_zero();
    // println!("mu_ik {:?} d_ij {:?}", mu_ik, d_ij/s_j);
    let mut x: f64 = 0.;

    loop {
        let calced_prob = prob_mu_ik_theta_i_x(x as f64, d_ij, mu_ik, t_ij, theta_i, s_j);
        if calced_prob > max_prob {
            max_prob = calced_prob;
        }
        result = result.ln_add_exp(calced_prob);

        if x > 200. && calced_prob - max_prob < LogProb(0.001_f64.ln()) {
            break;
        }
        x = x + 1.;
    }
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
    // (p1 - b + p2).exp() / (x + n)
    LogProb((p1 - b + p2) - (x + n).ln()) // TODO is this the correct form for returning LogProb?
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn test_neg_binom() {
        assert_relative_eq!(neg_binom(0., 10., 2.38).exp(), 0.2594752460369642);
        assert_relative_eq!(neg_binom(0., 30., 2.38).exp(), 0.16542351363026533);
        assert_relative_eq!(
            neg_binom(0., 30., 500.38).exp(),
            0.9809648435381609,
            epsilon = 1e-7
        );
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
