//! This implements formula 9 of the document and calculates the fold change / differential expression.
use std::collections::VecDeque;
use std::mem;
use std::path::Path;

use anyhow::Result;
use bio::stats::LogProb;
use rayon::prelude::*;

use crate::common::{Outdir, Pair};
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::sample_expression;


pub(crate) fn diff_exp(
    c: f64,
    preprocessing: &Path,
    group_path1: &Path,
    group_path2: &Path,
    out_dir: &Path,
) -> Result<()> {
    let out_dir = Outdir::create(out_dir)?;

    let in_dir1 = Outdir::open(&group_path1)?;
    let in_dir2 = Outdir::open(&group_path2)?;

    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190432).collect();

    let subsampled_ids = vec!["ENST00000671775.2", "ENST00000643797.1", "ENST00000496791.1", 
    "ENST00000651279.1", "ENST00000533357.5", "ENST00000538709.1", "ENST00000539934.5", 
    "ENST00000541259.1", "ENST00000542241.5", "ENST00000543262.5", "ENST00000549259.5", 
    "ENST00000551671.5", "ENST00000553379.6", "ENST00000553786.1", "ENST00000563608.2", "ENST00000566566.2"];

    // feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(_, feature_id)| -> Result<()> {
            // if subsampled_ids.contains(&feature_id.as_str()) {
            if feature_id.as_str() == "ERCC-00130" { 
            let prob_dist_i_k1: ProbDistribution2d = in_dir1.deserialize_value(feature_id)?;
            let prob_dist_i_k2: ProbDistribution2d = in_dir2.deserialize_value(feature_id)?;

            let mut diff_exp_distribution = ProbDistribution1d::new();

            let calc_prob = |f: f64| {
                    let calc_prob_fixed_theta = |theta|{
                    let density_x= |_, x: f64| {
                        prob_dist_i_k1.get(&[x, theta]) + prob_dist_i_k2.get(&[(f * (x + c) - c), theta])
                    };
                    // let mut points1 = prob_dist_i_k1.points.keys().map(|value| -c + (c + value.raw()) / f).collect::<Vec<_>>();
                    // let mut points2 = prob_dist_i_k2.points.keys().map(|value| value.raw()).collect::<Vec<_>>();
    
                    // points1.append(&mut points2);
                    // points1.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
                    // points1.dedup();              
    
                    // let probs_x = points1
                    //                 .windows(2)
                    //                 .map(|x| LogProb::ln_simpsons_integrate_exp(density_x, x[0], x[1], 5))
                    //                 .collect::<Vec<_>>();
                    // let prob_x= LogProb::ln_sum_exp(&probs_x);
                    let x_range = prob_dist_i_k1.get_range_per_theta(theta);
                    // println!("theta {:?}, x_range {:?}", theta, x_range);
                    let prob_x = LogProb::ln_simpsons_integrate_exp(density_x, 
                        x_range[0], 
                        x_range[1], 
                        301);
                    // println!("prob_x {:?}", prob_x.exp());
                    prob_x
                };
                
                let density_theta = |_, theta:f64|{
                    calc_prob_fixed_theta(theta) + LogProb(*prior.prob(theta))
                };

                let prob = LogProb::ln_simpsons_integrate_exp(
                        density_theta,
                        prior.min_value(),
                        prior.max_value(),
                        201,
                    );
                // let prob1 = LogProb::ln_simpsons_integrate_exp(
                //         density_theta,
                //         prior.min_value(),
                //         20.,
                //         151,
                // );
                // let prob2 = LogProb::ln_simpsons_integrate_exp(
                //     density_theta,
                //     20.,
                //     prior.max_value(),
                //     101,
                // );
                // let prob = prob1.ln_add_exp(prob2);
                
                println!("f {:?}, prob {:?}", f, prob.exp());
                // let prob = LogProb::ln_simpsons_integrate_exp(density, 0., 15., 11);
                prob
            };
            for x in 0..10{
                let f = x as f64 / 10. ;
                println!("f {:?}", f);
                diff_exp_distribution.insert(f as f64, calc_prob(f as f64));
            }
            for f in 1..11{
                println!("f {:?}", f);
                diff_exp_distribution.insert(f as f64, calc_prob(f as f64));

                // println!("f {:?}", f as f64 +0.5);
                // diff_exp_distribution.insert(f as f64 + 0.5, calc_prob(f as f64 +0.5));
            }

            // let mut start_points = vec![0.];
            // let mut cur_prob = LogProb::ln_zero();
            // diff_exp_distribution.insert(0., cur_prob);
            // let mut cur_max_prob_fold_change = 1. / max_prob_fold_change / 16.;
            // if cur_max_prob_fold_change > 0. {
            //     while cur_max_prob_fold_change < 32. {
            //         start_points.push(cur_max_prob_fold_change);
            //         cur_prob = calc_prob(cur_max_prob_fold_change);
            //         diff_exp_distribution.insert(cur_max_prob_fold_change, cur_prob);
            //         // println!("insert mu {:?}, prob {:?}", cur_maximum_likelihood_mean, f64::from(cur_prob.exp()));
            //         cur_max_prob_fold_change = cur_max_prob_fold_change * (2.0_f64).sqrt();
            //     }
            // }
            // if start_points.len() == 1 {
            //     start_points.push(5.);
            //     cur_prob = calc_prob(5.);
            //     diff_exp_distribution.insert(5., cur_prob);
            //     start_points.push(16.);
            //     cur_prob = calc_prob(16.);
            //     diff_exp_distribution.insert(16., cur_prob);
            //     // println!("insert mu {:?}, prob {:?}", 5000., f64::from(cur_prob.exp()));
            // }
            // start_points.push(32.);
            // cur_prob = calc_prob(32.);
            // diff_exp_distribution.insert(32., cur_prob);
            // println!("start points fs {:?}", start_points);


            // let mut queue = VecDeque::<Pair>::new();
            // start_points.windows(2).for_each(|w| {
            //     queue.push_back(Pair {
            //         left: w[0],
            //         right: w[1],
            //     })
            // });

            // while queue.len() > 0 {
            //     // println!("Len of queue {:?}", queue.len());
            //     // println!("#values inserted {:?}", diff_exp_distribution.len());
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
            //     let estimated_value = diff_exp_distribution.get(middle);
            //     let calculated_value = calc_prob(middle);
            //     // println!("middle {:?}, est: {:?}, calc: {:?}", middle, estimated_value.exp(), calculated_value.exp());
            //     diff_exp_distribution.insert(middle, calculated_value);
            //     if (estimated_value.exp() - calculated_value.exp()).abs() > 0.1 {
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
            
            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;
        }
            Ok(())
        })?;

    Ok(())
}
