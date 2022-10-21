//! This implements formula 3+4 of the document.
use std::fs;
use std::mem;
use std::path::Path;
// use std::process::exit;


use anyhow::Result;
use bio::stats::LogProb;
use getset::Getters;
// use itertools_num::linspace;
use rayon::prelude::*;
use rmp_serde::Deserializer;
use serde::Deserialize as SerdeDeserialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;

use crate::common::{window_x, Outdir}; // , Square, Point
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::ProbDistribution2d;

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
    let group_means = preprocessing.group_means();

    let s_j = preprocessing
        .scale_factors()
        .get(sample_id)
        .unwrap()
        .clone();

    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();
    
    let subsampled_ids = vec!["ENST00000671775.2", "ENST00000643797.1", "ENST00000496791.1", 
        "ENST00000651279.1", "ENST00000533357.5", "ENST00000538709.1", "ENST00000539934.5", 
        "ENST00000541259.1", "ENST00000542241.5", "ENST00000543262.5", "ENST00000549259.5", 
        "ENST00000551671.5", "ENST00000553379.6", "ENST00000553786.1", "ENST00000563608.2", "ENST00000566566.2"];
    let out_dir = Outdir::create(out_dir)?;
    feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            if subsampled_ids.contains(&feature_id.as_str()) {
            // if feature_id.as_str() == "ENST00000566566.2" {   

            let d_ij = mean_disp_estimates.means()[*i]; //TODO Do we need group mean mu_ik instead of sample mean
            let mut d_ik = group_means[*i];
            if d_ik < 1e-8 { //TODO SEnsible? Filter at different position better?
                d_ik = 0.;
            }
            // METHOD: If the per-sample dispersion is unknown, fall back to a mean interpolated from the other samples.
            let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                // TODO log message
                return Ok(());
            };
            // println!("\n\n -------------------------- feature {:?}", i);
            // println!("d_ij {:?}, d_ik {:?}, t_ij {:?}", d_ij, d_ik, t_ij);

            let mut likelihoods = ProbDistribution2d::new();
            // println!("--------------------------3");
            let mut start_points_mu_ik = vec![0.];
            let mut cur_d_ik = d_ik / 16.;
            if (cur_d_ik > 0.) {
                while cur_d_ik < 10000. {
                    start_points_mu_ik.push(cur_d_ik);
                    cur_d_ik = cur_d_ik * 2.;
                }
            }
            if start_points_mu_ik.len() == 1 {
                start_points_mu_ik.push(5000.);
            }
            start_points_mu_ik.push(10001.);


            // let mut start_points_mu_ik = vec![0.];
            // if d_ik != 0. {
            //     let mut cur_d_ik = d_ik / 10.;
            //     for i in 1..4  {
            //         if cur_d_ik >= 10000. { break; }
            //         start_points_mu_ik.push(cur_d_ik);
            //         cur_d_ik = cur_d_ik * 10.;
            //     }
            //     if start_points_mu_ik.len() == 1 {
            //         start_points_mu_ik.push(5000.);
            //     }
            //     start_points_mu_ik.push(10000.);
            // }
            // println!("start_points_mu_ik {:?}",start_points_mu_ik);
            let mut start_points_theta_i = Vec::<f64>::new();
            start_points_theta_i.extend(prior.left_window());
            start_points_theta_i.extend(prior.right_window());
            start_points_theta_i.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let calc_prob = |m, t| {
                // println!("XXX mu {:?}, theta {:?}", m, t);
                likelihood_mu_ik_theta_i(
                    d_ij,
                    m,  // mu_ik
                    t_ij,
                    t,  // theta_i
                    s_j,
                    epsilon,
                )
            };
            likelihoods.insert_grid(start_points_mu_ik, start_points_theta_i, calc_prob);

            
            // let theta_i = 1000.;
            // for mu_ik in [1., 100., 500., 1000., 2500., 5000.] {
            //     let mut result = LogProb::ln_zero();
            //     println!("mu {:?}, theta {:?}, s_j {:?}",mu_ik, theta_i, s_j);
            //     for x in [0.,1.,10., 50.,100., 500., 1000., 5000., 10000., 50000.]{
            //         let nb1 = neg_binom(d_ij, x, t_ij);
            //         let nb2 = neg_binom(x, mu_ik * s_j, theta_i);
            //         let sum = nb1 + nb2;
            //         // result = result.ln_add_exp(sum);
            //         println!("{:?} & {:?} & {:?} & {:?} \\\\", x, f64::from(nb1).exp(), f64::from(nb2).exp(), f64::from(sum).exp() );
            //     }
            //     for x in 0..50001 {
            //         let nb1 = neg_binom(d_ij, x as f64, t_ij);
            //         let nb2 = neg_binom(x as f64, mu_ik * s_j, theta_i);
            //         let sum = nb1 + nb2;
            //         result = result.ln_add_exp(sum);
            //     }
            //     println!("result {:?}", f64::from(result).exp());
            // }
            

            out_dir.serialize_value(feature_id, likelihoods)?;
            }
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
    // if d_ij == 0. {
    //     // println!("############ prob_mu_ik_theta_i_x x {:?}, d_ij {:?}, mu_ik {:?}, t_ij {:?}, theta_i {:?}, s_j {:?}", x, d_ij, mu_ik, t_ij, theta_i, s_j);
    //     println!("x {:?}, nb(d_ij, x, t_ij) {:?}, nb(x, mu_ik * s_j, theta_i) {:?}, res {:?}",x, neg_binom(d_ij, x, t_ij), neg_binom(x, mu_ik * s_j, theta_i),neg_binom(d_ij, x, t_ij) + neg_binom(x, mu_ik * s_j, theta_i));
    // }
    neg_binom(d_ij, x, t_ij) + neg_binom(x, mu_ik * s_j, theta_i)
}

/// Inner of equation 3/4 in the document.
fn likelihood_mu_ik_theta_i(
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    s_j: f64,
    _: LogProb,   // epsilon
) -> LogProb {
    if d_ij != 0. && mu_ik == 0. {
        return LogProb::ln_zero();
    }
    if mu_ik == 3208600050.4032364 {
        println!("############ likelihood_mu_ik_theta_i");
        println!("mu_ik {:?}, theta_i {:?}, d_ij {:?}, t_ij {:?}, s_j {:?}", mu_ik, theta_i, d_ij, t_ij, s_j);
    }
    let mut max_prob = LogProb::ln_zero();
    let mut result = LogProb::ln_zero();
    let mut x : u64 = 0;
    while true{
        let calced_prob = prob_mu_ik_theta_i_x(x as f64, d_ij, mu_ik, t_ij, theta_i, s_j);
        if calced_prob > max_prob{
            max_prob = calced_prob;
        }       
        result = result.ln_add_exp(calced_prob);
        if x > 500 && calced_prob - max_prob < LogProb(0.01_f64.ln()){
            break;
        }
        x = x+1;
    }
    // println!("mu_ik {:?}, theta_i {:?}, x {:?}", mu_ik, theta_i, x);
    result
}

fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);

    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    // (p1 - b + p2).exp() / (x + n)
    if x == 0. {
        // println!("x {:?}, mu {:?}, theta {:?}, n {:?}, p {:?}, p1 {:?}, p2 {:?}, b {:?}", x, mu, theta, n, p ,p1, p2, b);
    }
    LogProb((p1 - b + p2) - (x + n).ln()) // TODO is this the correct form for returning LogProb?
}

#[cfg(test)]
mod tests{
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn test_neg_binom(){
        assert_relative_eq!(neg_binom(0., 10., 2.38).exp(), 0.2594752460369642);
        assert_relative_eq!(neg_binom(0., 30., 2.38).exp(), 0.16542351363026533);
        assert_relative_eq!(neg_binom(0., 30., 500.38).exp(), 0.9809648435381609, epsilon=1e-7);
        
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
