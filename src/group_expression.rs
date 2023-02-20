//! This implements formula 5,6,7 of the document.


use std::path::{Path, PathBuf};
use std::collections::VecDeque;
use std::fs::File;

use anyhow::Result;
use bio::stats::LogProb;
use rayon::prelude::*;
use csv;
use itertools_num::linspace;

use crate::common::{QueryPoints, Outdir, Pair, difference_to_big};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::sample_expression::SampleInfo;



pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    c: f64,
    out_dir: &Path,
) -> Result<()> {
    let query_points = QueryPoints::new(c)?;
    let mu_ik_points = query_points.all_mu_ik();
    let start_points_theta_i = query_points.thetas();
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190400).collect(); //190432

    let file = File::open("/vol/nano/bayesian-diff-exp-analysis/texp-evaluation/estimated_dispersion.csv")?;
    let mut rdr = csv::Reader::from_reader(file);
    let mut thetas = Vec::<f64>::new();
    for result in rdr.records() {
        let record = result?;
        let dispersion: f64 = record[1].parse().unwrap();
        thetas.push(dispersion);
    }


    let subsampled_ids = vec!["ERCC-00130","ERCC-00004", "ERCC-00136", "ERCC-00096", "ERCC-00171", "ERCC-00009",
    "ERCC-00074", "ERCC-00113", "ERCC-00145", "ERCC-00002", "ERCC-00046", "ERCC-00003"];

    let out_dir = Outdir::create(out_dir)?;
    // feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            // if subsampled_ids.contains(&feature_id.as_str()) {
            // if feature_id.as_str() == "ERCC-00130" {    
            
            println!("--------------feature {:?} {:?}", i, feature_id);
            let maximum_likelihood_means: Vec<f64> = sample_expression_paths
                .iter()
                .map(|sample_expression_path| {
                    let sample_info: SampleInfo =
                        Outdir::open(sample_expression_path)?.deserialize_value("info")?;
                    Ok(preprocessing
                        .mean_disp_estimates()
                        .get(sample_info.sample_id())
                        .ok_or(Error::UnknownSampleId {
                            sample_id: sample_info.sample_id().to_owned(),
                        })?
                        .means()[*i])
                })
                .collect::<Result<Vec<_>>>()?;
            // println!("1");
            let sample_expression_likelihoods = sample_expression_paths
                .iter()
                .map(|sample_expression_path| {
                    let dir = Outdir::open(sample_expression_path)?;
                    let feature_id_with_mpk = format!("{}{}", feature_id, ".mpk");
                    let mut fullpath = format!(
                        "{}{}",
                        sample_expression_path.to_str().unwrap(),
                        feature_id_with_mpk
                    );
                    if sample_expression_path
                        .to_str()
                        .unwrap()
                        .chars()
                        .last()
                        .unwrap()
                        != '/'
                    {
                        fullpath = format!(
                            "{}/{}",
                            sample_expression_path.to_str().unwrap(),
                            feature_id_with_mpk
                        );
                    }

                    if Path::new(&fullpath).exists() {
                        let likelihoods: ProbDistribution2d = dir.deserialize_value(feature_id)?;
                        Ok(likelihoods)
                    } else {
                        Ok(ProbDistribution2d::na())
                    }
                })
                .collect::<Result<Vec<_>>>()?;
                if sample_expression_likelihoods.iter().all(|x| x.is_na()){
                    out_dir.serialize_value(feature_id, ProbDistribution2d::na())?;
                    return Ok(());
                }
            // println!("2");
            let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
                / maximum_likelihood_means.len() as f64;

            // let mut prob_dist = ProbDistribution1d::new();
            let mut prob_dist = ProbDistribution2d::new();
            // println!("3");
            let calc_prob = |mu_ik : f64, theta_i: f64| {
                // println!("mu_ik {:?}", mu_ik);
                if mu_ik == 0. {
                    return LogProb::ln_zero();
                }
                let mut prob = sample_expression_likelihoods
                    .iter()
                    .map(|sample_expression_likelihood| {
                        sample_expression_likelihood
                        .get(&[mu_ik, theta_i])
                    })
                    .sum::<LogProb>();  //Formula 5
                    // +LogProb(*prior.prob(theta_i));
                // prob = LogProb(f64::from(prob) * 8.);

                // Result of formula 7.
                // let prob= LogProb::ln_simpsons_integrate_exp(
                //     density,
                //     prior.min_value(),
                //     prior.max_value(),
                //     451,
                // );
                // let prob = density(0., prior.mean());
                // println!("mu_ik {:?}, prob {:?}", mu_ik, prob);
                prob

            };
            // let mut start_points_mu_ik = vec![0.];
            // // let mut cur_prob = calc_prob(0.);
            // // prob_dist.insert(0., cur_prob);            
            // // // println!("insert mu {:?}, prob {:?}", 0., f64::from(cur_prob.exp()));
            // let mut cur_maximum_likelihood_mean = maximum_likelihood_mean / 16.;
            // if (cur_maximum_likelihood_mean > 0.) {
            //     while cur_maximum_likelihood_mean < 10000. {
            //         start_points_mu_ik.push(cur_maximum_likelihood_mean);
            // //         cur_prob = calc_prob(cur_maximum_likelihood_mean);
            // //         prob_dist.insert(cur_maximum_likelihood_mean, cur_prob);
            // //         // println!("insert mu {:?}, prob {:?}", cur_maximum_likelihood_mean, f64::from(cur_prob.exp()));
            //         cur_maximum_likelihood_mean = cur_maximum_likelihood_mean * 2.;
            //     }
            // }
            // if start_points_mu_ik.len() == 1 {
            //     start_points_mu_ik.push(2500.);
            // //     cur_prob = calc_prob(5000.);
            // //     prob_dist.insert(5000., cur_prob);
            // //     // println!("insert mu {:?}, prob {:?}", 5000., f64::from(cur_prob.exp()));
            // }
            // start_points_mu_ik.push(5000.);
            // // cur_prob = calc_prob(10000.);
            // // prob_dist.insert(10000., cur_prob);
            // // // println!("insert mu {:?}, prob {:?}", 10000., f64::from(cur_prob.exp()));
            let mut start_points_mu_ik = mu_ik_points.clone();
            
            // let theta = thetas[i-190432];

            // let mut start_points_theta_i = Vec::<f64>::new();
            // start_points_theta_i.extend(prior.left_window());
            // start_points_theta_i.push(theta);
            // start_points_theta_i.extend(prior.right_window());
            // start_points_theta_i.sort_by(|a, b| a.partial_cmp(b).unwrap());
           



            // let mut queue = VecDeque::<Pair>::new();
            // start_points.windows(2).for_each(|w| {
            //     queue.push_back(Pair {
            //         left: w[0],
            //         right: w[1],
            //     })
            // });
            // while queue.len() > 0 {
            //     // println!("--------------------------");
            //     let pair = queue.pop_front().unwrap();
            //     let left = pair.left;
            //     let right = pair.right;
            //     if left == 0. && right == 0. {
            //         continue;
            //     }
            //     let mut middle = 0.;
            //     if left.is_finite() && right.is_finite() {
            //         middle = left / 2. + right / 2.;
            //     } else if left.is_finite() {
            //         middle = 10. * left;
            //     }
            //     if middle.is_infinite() {
            //         continue;
            //     }
            //     let estimated_value = prob_dist.get(middle);
            //     let calculated_value = calc_prob(middle);
            //     // println!("middle {:?},est {:?}, calc {:?}, len {:?}", middle, estimated_value, calculated_value, queue.len());
            //     prob_dist.insert(middle, calculated_value);
            //     // println!("insert mu {:?}, prob {:?}", middle, f64::from(calculated_value.exp()));
            //     // println!("left {:?}, middle {:?}, right {:?}", left, middle, right);
            //     // println!("est {:?}, calc {:?}", estimated_value, calculated_value);
            //     let prob_dist_max = prob_dist.get_max_prob();
            //     if middle - left < 1e-3 || prob_dist.len() > 10000{
            //         continue;
            //     }
            //     if difference_to_big(estimated_value, calculated_value, prob_dist_max) {
            //         queue.push_back(Pair {
            //             left: left,
            //             right: middle,
            //         });
            //         queue.push_back(Pair {
            //             left: middle,
            //             right: right,
            //         });
            //     }
            // }

            // println!("4");
            prob_dist.insert_grid(start_points_mu_ik, start_points_theta_i.clone(), calc_prob);
            // println!("5");
            // let norm_factor = prob_dist.normalize(); // remove factor c_ik
            out_dir.serialize_value(feature_id, prob_dist)?;
            // }
            Ok(())
        })?;

    Ok(())
}
