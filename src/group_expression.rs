//! This implements formula 5,6,7 of the document.


use std::path::{Path, PathBuf};
use std::collections::VecDeque;

use anyhow::Result;
use bio::stats::LogProb;
use rayon::prelude::*;

use crate::common::{Outdir, Pair, difference_to_big};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::sample_expression::SampleInfo;



pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    out_dir: &Path,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
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
            
            // println!("FEATURE --------------------------");
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
            let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
                / maximum_likelihood_means.len() as f64;

            let mut prob_dist = ProbDistribution1d::new();

            let calc_prob = |mu_ik : f64| {
                // println!("mu_ik {:?}", mu_ik);
                if mu_ik == 0. {
                    return LogProb::ln_zero();
                }
                let density = |_, theta_i| {
                    let d = sample_expression_likelihoods
                        .iter()
                        .map(|sample_expression_likelihood| {
                            sample_expression_likelihood
                                .get(&[mu_ik, theta_i])
                        })
                        .sum::<LogProb>() + //Formula 5
                        LogProb(*prior.prob(theta_i) * 2.0); // square of Pr(theta_i), formula 8
                    d
                };

                // Result of formula 7.
                let prob= LogProb::ln_simpsons_integrate_exp(
                    density,
                    prior.min_value(),
                    prior.max_value(),
                    451,
                );
                prob

            };
            // println!("\n i {:?}, maximum_likelihood_mean {:?}", i, maximum_likelihood_mean);
            let mut start_points = vec![0.];
            let mut cur_prob = calc_prob(0.);
            prob_dist.insert(0., cur_prob);            
            // println!("insert mu {:?}, prob {:?}", 0., f64::from(cur_prob.exp()));
            let mut cur_maximum_likelihood_mean = maximum_likelihood_mean / 16.;
            if (cur_maximum_likelihood_mean > 0.) {
                while cur_maximum_likelihood_mean < 10000. {
                    start_points.push(cur_maximum_likelihood_mean);
                    cur_prob = calc_prob(cur_maximum_likelihood_mean);
                    prob_dist.insert(cur_maximum_likelihood_mean, cur_prob);
                    // println!("insert mu {:?}, prob {:?}", cur_maximum_likelihood_mean, f64::from(cur_prob.exp()));
                    cur_maximum_likelihood_mean = cur_maximum_likelihood_mean * 2.;
                }
            }
            if start_points.len() == 1 {
                start_points.push(5000.);
                cur_prob = calc_prob(5000.);
                prob_dist.insert(5000., cur_prob);
                // println!("insert mu {:?}, prob {:?}", 5000., f64::from(cur_prob.exp()));
            }
            start_points.push(10000.);
            cur_prob = calc_prob(10000.);
            prob_dist.insert(10000., cur_prob);
            // println!("insert mu {:?}, prob {:?}", 10000., f64::from(cur_prob.exp()));
            
            let mut queue = VecDeque::<Pair>::new();
            start_points.windows(2).for_each(|w| {
                queue.push_back(Pair {
                    left: w[0],
                    right: w[1],
                })
            });
            while queue.len() > 0 {
                // println!("--------------------------");
                let pair = queue.pop_front().unwrap();
                let left = pair.left;
                let right = pair.right;
                if left == 0. && right == 0. {
                    continue;
                }
                let mut middle = 0.;
                if left.is_finite() && right.is_finite() {
                    middle = left / 2. + right / 2.;
                } else if left.is_finite() {
                    middle = 10. * left;
                }
                if middle.is_infinite() {
                    continue;
                }
                let estimated_value = prob_dist.get(middle);
                let calculated_value = calc_prob(middle);
                // println!("middle {:?},est {:?}, calc {:?}, len {:?}", middle, estimated_value, calculated_value, queue.len());
                prob_dist.insert(middle, calculated_value);
                // println!("insert mu {:?}, prob {:?}", middle, f64::from(calculated_value.exp()));
                // println!("left {:?}, middle {:?}, right {:?}", left, middle, right);
                // println!("est {:?}, calc {:?}", estimated_value, calculated_value);
                let prob_dist_max = prob_dist.get_max_prob();
                if middle - left < 1e-3 || prob_dist.len() > 10000{
                    continue;
                }
                if difference_to_big(estimated_value, calculated_value, prob_dist_max) {
                    queue.push_back(Pair {
                        left: left,
                        right: middle,
                    });
                    queue.push_back(Pair {
                        left: middle,
                        right: right,
                    });
                }
            }

            let norm_factor = prob_dist.normalize(); // remove factor c_ik
            out_dir.serialize_value(feature_id, prob_dist)?;
            }
            Ok(())
        })?;

    Ok(())
}
