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
    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();
    feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(_, feature_id)| -> Result<()> {
            // if feature_id.as_str() == "ENST00000641952.1" {
            let mut prob_dist_i_k1: ProbDistribution1d = in_dir1.deserialize_value(feature_id)?;
            let mut prob_dist_i_k2: ProbDistribution1d = in_dir2.deserialize_value(feature_id)?;

            let mut max_prob_position1 = prob_dist_i_k1.get_max_prob_position();
            if max_prob_position1 == 0. {
                max_prob_position1 = 1.;
            }
            let mut max_prob_position2 = prob_dist_i_k2.get_max_prob_position();
            if max_prob_position2 == 0. {
                max_prob_position2 = 1.;
            }; // Avoid division by zero
               // let max_prob_fold_change = (prob_dist_i_k1.get(max_prob_position1) - prob_dist_i_k2.get(max_prob_position2)).exp();
            

            if max_prob_position1 > max_prob_position2 {
                // Ensure that maximum of prob_dist_k1 is left of prob_dist_k2
                mem::swap(&mut prob_dist_i_k1, &mut prob_dist_i_k2);
                mem::swap(&mut max_prob_position1, &mut max_prob_position2);
            }
            let max_prob_fold_change = max_prob_position1 / max_prob_position2;

            let mut diff_exp_distribution = ProbDistribution1d::new();

            if max_prob_fold_change < 0.00001 {
                out_dir.serialize_value(feature_id, diff_exp_distribution)?;
                return Ok(());
            }

            let calc_prob = |f: f64| {
                // println!("f {:?}", f);
                // let f_max_i_k2 = (-f * c + c + max_prob_position1) / f;
                let density = |_, x: f64| {
                    // println!("prob 1 {:?} von x {:?}", prob_dist_i_k1.get(x), x);
                    // println!("prob 2 {:?} von x * f {:?}", prob_dist_i_k2.get(f * (x + c) - c), f*x);
                    prob_dist_i_k1.get(x) + (prob_dist_i_k2.get(f * (x + c) - c))
                };
                let mut points1 = prob_dist_i_k1.points.keys().map(|value| value.raw()).collect::<Vec<_>>();
                let mut points2 = prob_dist_i_k2.points.keys().map(|value| -c + (c + value.raw()) / f).collect::<Vec<_>>();
                points1.append(&mut points2);
                points1.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
                points1.dedup();               
                // println!("num points {:?}", points1.len());
                let probs = points1
                                .windows(2)
                                .map(|x| LogProb::ln_simpsons_integrate_exp(density, x[0], x[1], 3))
                                .collect::<Vec<_>>();
                let prob = LogProb::ln_sum_exp(&probs);

                
                // let prob = LogProb::ln_simpsons_integrate_exp(density, 0., 15., 11);
                prob
            };

            // println!("max_prob_fold_change {:?}", max_prob_fold_change);
            let start_points = [
                0.,
                max_prob_fold_change / 10.,
                max_prob_fold_change,
                max_prob_fold_change * 10.,
                f64::INFINITY,
            ];
            // println!("Len of start points mu_ik {:?}", start_points.len());
            diff_exp_distribution.insert(
                max_prob_fold_change / 10.,
                calc_prob(max_prob_fold_change / 10.),
            );
            diff_exp_distribution.insert(max_prob_fold_change, calc_prob(max_prob_fold_change));
            diff_exp_distribution.insert(
                max_prob_fold_change * 10.,
                calc_prob(max_prob_fold_change * 10.),
            );
            let mut queue = VecDeque::<Pair>::new();
            start_points.windows(2).for_each(|w| {
                queue.push_back(Pair {
                    left: w[0],
                    right: w[1],
                })
            });

            while queue.len() > 0 {
                // println!("Len of queue {:?}", queue.len());
                // println!("#values inserted {:?}", diff_exp_distribution.len());
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
                let estimated_value = diff_exp_distribution.get(middle);
                // println!("est: {:?}", estimated_value);
                let calculated_value = calc_prob(middle);
                // println!("calc: {:?}", calculated_value);
                diff_exp_distribution.insert(middle, calculated_value);
                if (estimated_value.exp() - calculated_value.exp()).abs() > 0.01 {
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

            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;
            Ok(())
        })?;

    Ok(())
}
