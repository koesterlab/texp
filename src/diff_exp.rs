//! This implements formula 9 of the document and calculates the fold change / differential expression.

use std::mem;
use std::path::Path;
use std::collections::VecDeque;
use std::collections::LinkedList;

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use noisy_float::types::N64;
use rayon::prelude::*;

use crate::common::{window_f, Log2FoldChange, Mean, Outdir};
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::preprocess::Preprocessing;


pub(crate) struct Pair {
    left: f64,
    right: f64,
}

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
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            let mut prob_dist_i_k1: ProbDistribution1d = in_dir1.deserialize_value(feature_id)?;
            let mut prob_dist_i_k2: ProbDistribution1d = in_dir2.deserialize_value(feature_id)?;

            let mut max_prob_value1 = prob_dist_i_k1.get_max_prob_value();
            let mut max_prob_value2 = prob_dist_i_k2.get_max_prob_value();
            let mut max_prob_value2 = if max_prob_value2 != 0. {
                max_prob_value2
            } else {
                1.
            }; // Avoid division by zero
            // let max_prob_fold_change = (prob_dist_i_k1.get(max_prob_value1) - prob_dist_i_k2.get(max_prob_value2)).exp();
            let max_prob_fold_change = max_prob_value1 / max_prob_value2;

            if max_prob_value1 < max_prob_value2 {
                // Ensure that maximum of prob_dist_k1 is left of prob_dist_k2
                mem::swap(&mut prob_dist_i_k1, &mut prob_dist_i_k2);
                mem::swap(&mut max_prob_value1, &mut max_prob_value2);
            }

            // Step 1: use window() to determine range around max_prob_fold_change
            let (left_window, right_window) = window_f(max_prob_fold_change);

            let mut diff_exp_distribution = ProbDistribution1d::new();

            let mut calc_prob = |f: f64| {
                let f_max_i_k2 = (-f * c + c + max_prob_value1) / f;
                // for x, iterate over values in prob_dist_i_k1, the value for looking up in prob_dist_i_k2 is f * (x + c) - c
                // let prob = LogProb::ln_trapezoidal_integrate_grid_exp(
                //     |i, &value| {
                //         prob_dist_i_k1.get(&value)
                //             + prob_dist_i_k2.get(&[(N32::new(f as f32) * (value[0] + N32::new(c))- N32::new(c))])
                //     },
                // &prob_dist_i_k1
                //     .points
                //     .keys()
                //     .map(|value| **value)
                //     .collect::<Vec<_>>(),
                // );

                let density = |i: usize, x: f64| {
                    prob_dist_i_k1
                        .get(x)
                        + (prob_dist_i_k2.get((f * (x + c) - c)))
                };

                let mut prob = LogProb::ln_simpsons_integrate_exp(density, 0., 15., 11);

                // let prob = LogProb::ln_one();
                // diff_exp_distribution.insert(f, prob);

                prob
            };
            println!("max_prob_fold_change {:?}", max_prob_fold_change);
            let start_points = [0., max_prob_fold_change/10., max_prob_fold_change, max_prob_fold_change*10., f64::INFINITY ];
            diff_exp_distribution.insert(max_prob_fold_change/10., calc_prob( max_prob_fold_change/10.));
            diff_exp_distribution.insert(max_prob_fold_change, calc_prob( max_prob_fold_change));
            diff_exp_distribution.insert(max_prob_fold_change*10.,calc_prob( max_prob_fold_change*10.));
            let mut queue = VecDeque::<Pair>::new(); 
            start_points.windows(2).for_each(|w| queue.push_back(Pair { left: w[0], right: w[1] }));
            while queue.len() > 0 {
                let pair = queue.pop_front().unwrap();
                let left = pair.left;
                let right = pair.right;
                if left == 0. && right == 0.{
                    continue
                }
                let mut middle = 0.;
                if left.is_finite() && right.is_finite(){
                    middle = left/2. + right/2.;
                }
                else if left.is_finite() {
                    middle = 10. * left;                    
                }
                if middle.is_infinite(){
                    continue
                }
                // println!("middle {:?}", middle);
                let estimated_value = diff_exp_distribution.get(middle);
                let calculated_value = calc_prob(middle);
                diff_exp_distribution.insert(middle, calculated_value);
                if (estimated_value.exp() - calculated_value.exp()).abs() > 0.01 {
                    queue.push_back(Pair { left: left, right: middle });
                    queue.push_back(Pair { left: middle, right: right });
                    // if left > 1e10 {
                    //     println!("left {:?}, right {:?}, middle {:?}, estimated_value {:?}, calculated_value {:?}", left, right, middle, estimated_value, calculated_value);
                    //     println!("est-calc {:?}, .abs {:?}",(estimated_value - calculated_value), (estimated_value - calculated_value).abs() )
                    // }
                }
                
            }
            // // Step 2: For each fold change in determined range, calculate formula 9 and put in ProbabilityDistribution
            // for f in left_window {
            //     calc_prob(f);
            // }

            // for f in right_window {
            //     calc_prob(f);
            // }

            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;

            Ok(())
        })?;

    Ok(())
}
