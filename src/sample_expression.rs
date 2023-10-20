//! This implements formula 3+4 of the document.
use std::fs;
use std::mem;
use std::path::Path;
// use std::process::exit;
// use std::fs::File;

// use csv;
use anyhow::Result;
use bio::stats::LogProb;
// use statrs::distribution::{Discrete,Poisson};
// use bio::stats::probs::adaptive_integration::ln_integrate_exp;
// use ordered_float::NotNan;
use getset::Getters;
use rayon::prelude::*;
use rmp_serde::Deserializer;
use serde::Deserialize as SerdeDeserialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;
// use itertools_num::linspace;

use crate::common::Outdir; 
// , Square, Point
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;

pub(crate) fn sample_expression(
    preprocessing: &Path,
    sample_id: &str,
    epsilon: LogProb,
    c : f64,
    out_dir_path: &Path,
) -> Result<()> {
    
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let sample_ids = preprocessing.scale_factors().keys().cloned().collect::<Vec<_>>();
    let prior = preprocessing.prior()?;
    let mean_disp_estimates =
        preprocessing
            .mean_disp_estimates()
            .get(sample_id)
            .ok_or(Error::UnknownSampleId {
                sample_id: sample_id.to_owned(),
            })?;

    let mut s_j = preprocessing
        .scale_factors()
        .get(sample_id)
        .unwrap()
        .clone();
    println!("s_j {:?}", s_j);
    // let temp: Vec<_>= mean_disp_estimates.means().iter().enumerate().skip(190432).collect();
    // println!("mean_disp_estimates {:?}", temp);

    let query_points = query_points::calc_query_points(c, preprocessing.mean_disp_estimates().clone(), sample_ids, preprocessing.feature_ids().clone());
    // let query_points = preprocessing.query_points();
    // let mu_ik_points = Vec::<f64>::new();//query_points.all_mu_ik();
    // let start_points_theta_i = Vec::<f64>::new(); //query_points.thetas();

    // let file = File::open("/vol/nano/bayesian-diff-exp-analysis/texp-evaluation/estimated_dispersion.csv")?;
    // let mut rdr = csv::Reader::from_reader(file);
    // let mut thetas = Vec::<f64>::new();
    // for result in rdr.records() {
    //     let record = result?;
    //     let dispersion: f64 = record[1].parse().unwrap();
    //     thetas.push(dispersion);
    // }



    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190432).collect();
    // println!("{:?} features", feature_ids.len());
    // println!("{:?}", feature_ids);
    // let subsampled_ids = vec!["ERCC-00147","ERCC-00117", "ERCC-00016", "ERCC-00086", "ERCC-0061", "ERCC-00048"];
    let subsampled_ids = vec!["ERCC-00130","ERCC-00004", "ERCC-00136", "ERCC-00096", "ERCC-00074", "ERCC-00002"];
    // ,    "ERCC-00171", "ERCC-00113", "ERCC-00145", "ERCC-00002", "ERCC-00046", "ERCC-00003", "ERCC-00113"];

    let out_dir = Outdir::create(out_dir_path)?;
    // feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            
            // if subsampled_ids.contains(&feature_id.as_str()) {
            // if feature_id.as_str() == "ERCC-00108" {   
            // prit start time of feature
            let time1= std::time::SystemTime::now();
            // println!("feature {:?} {:?} started at {:?}", i, feature_id, time1 );

            // println!("\n--------------feature {:?} {:?}", i, feature_id);
            // let mut output = out_dir_path.to_path_buf();
            // output.push(feature_id);
            // output.set_extension("csv");
            // let mut wtr = csv::Writer::from_path(output)?;
            // wtr.serialize(("mu_ik", "probability")).unwrap();

            let d_ij = mean_disp_estimates.means()[*i]; //TODO Do we need group mean mu_ik instead of sample mean
            // println!("{:?} d_ij {:?}", feature_id, d_ij / s_j);
            // METHOD: If the per-sample dispersion is unknown, fall back to a mean interpolated from the other samples.
            let mut t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                println!("skipped {:?}", feature_id);
                // TODO log message
                return Ok(());
            };
            // println!("t_ij {:?}", t_ij);
            
            
            // println!("d_ij {:?}, t_ij {:?}", d_ij, t_ij);

            let mut likelihoods = ProbDistribution2d::new();
            // println!("--------------------------3");

            

            // let mut start_points_mu_ik = vec![0.];
            // let mut cur_d_ik = d_ik / 16.;
            // if (cur_d_ik > 0.) {
            //     while cur_d_ik < 10000. {
            //         start_points_mu_ik.push(cur_d_ik);
            //         cur_d_ik = cur_d_ik * 2.;
            //     }
            // }
            // if start_points_mu_ik.len() == 1 {
            //     start_points_mu_ik.push(2500.);
            // }
            // start_points_mu_ik.push(5001.);

            // let theta = thetas[i-190432];
            // let mut start_points_theta_i: Vec<f64> = Vec::<f64>::new();

            // let mut start_points_theta_i: Vec<f64> = linspace(0.1, 1., 10).collect();
            // start_points_theta_i.extend( linspace(1.5, 10., 18));
            // start_points_theta_i.extend( linspace(11., 165., 155));
            // // println!("start_points_theta_i {:?}", start_points_theta_i);

            let calc_prob = |m, t| {
                // println!("mu {:?}, theta {:?}", m, t);
                let prob = likelihood_mu_ik_theta_i(
                    d_ij,
                    m,  // mu_ik
                    t_ij,
                    t,  // theta_i
                    s_j,
                    epsilon,
                );
                // if t == 0.01 {
                //     wtr.serialize((m, prob.exp())).unwrap();
                //     // if m < 5. {
                //     //     println!("mu {:?}, prob {:?}", m, prob.exp());
                //     // }
                // }
                prob
            };
            let mu_ik_points = query_points.get(&feature_id.to_string()).unwrap().all_mu_ik();
            let start_points_theta_i = query_points.get(&feature_id.to_string()).unwrap().thetas();
            // println!("{:?} #mu_ik {:?}, #theta_i {:?}",feature_id,  mu_ik_points.len(), start_points_theta_i.len());
            // println!("{:?} mu_ik_points {:?} {:?}",feature_id, mu_ik_points[0], mu_ik_points[mu_ik_points.len()-1]);
            // println!("start_points_theta_i {:?}", start_points_theta_i);
            likelihoods.insert_grid(d_ij/s_j, mu_ik_points.clone(), start_points_theta_i.clone(), calc_prob);

            out_dir.serialize_value(feature_id, likelihoods)?;
            let time2= std::time::SystemTime::now();
            println!("sample {} feature {:?} {:?} finished in duration {:?}",sample_id, i, feature_id, time2.duration_since(time1).unwrap());
            // }
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
    // println!("x {:?}, d_ij {:?}, mu_ik {:?}, t_ij {:?}, theta_i {:?}, s_j {:?}", x, d_ij, mu_ik, t_ij, theta_i, s_j);
    //METHOD /s_j?? mu_ik*s_j würde einluss von t_ij ändern, weil mu geändert wird.
    let left = neg_binom(d_ij/s_j, x, t_ij);
    let right = neg_binom(x, mu_ik, theta_i);
    // let right =  LogProb::from(Poisson::new(mu_ik).unwrap().ln_pmf(x as u64));
    let result = left + right; 
    // if theta_i == 0.01 && mu_ik < 5.{
    //     println!("x {:?}, mu_ik {:?}, left {:?}, right {:?}, result {:?}", x, mu_ik, left.exp(), right.exp(), result.exp());
    // }
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
    _: LogProb,   // epsilon
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
    let mut x : f64 = 0.;
    
    loop {
        let calced_prob = prob_mu_ik_theta_i_x(x as f64, d_ij, mu_ik, t_ij, theta_i, s_j);
        if calced_prob > max_prob{
            max_prob = calced_prob;
        }       
        result = result.ln_add_exp(calced_prob);

        if (x > 200. && calced_prob - max_prob < LogProb(0.001_f64.ln())) {
            break;
        }
        // if x % 1000. == 0. {
        //     println!("x {:?}, calced_prob {:?}, max_prob {:?} diff {:?}", x, calced_prob.exp(), max_prob.exp(), (calced_prob - max_prob).exp());
        // }
        x = x + 1.;
    }
    // println!("mu_ik {:?} x {:?}", mu_ik, x); 
    // if theta_i < 0.3 {
    //     println!("mu_ik {:?}, theta_i {:?}, x {:?}, result {:?}, max_prob {:?}", mu_ik, theta_i, x, result, max_prob);
    // }
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
