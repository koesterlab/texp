//! This implements formula 9 of the document and calculates the fold change / differential expression.

use std::path::Path;
use std::mem;

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rayon::prelude::*;

use crate::common::{window_f, Log2FoldChange, Mean, Outdir, ProbDistribution1d};
use crate::preprocess::Preprocessing;

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
            let max_prob_fold_change = max_prob_value1 / max_prob_value2;

            if max_prob_value1 < max_prob_value2 { // Ensure that maximum of prob_dist_k1 is left of prob_dist_k2
                mem::swap(&mut prob_dist_i_k1, &mut prob_dist_i_k2);
                mem::swap(&mut max_prob_value1, &mut max_prob_value2);
            }

            // Step 1: use window() to determine range around max_prob_fold_change
            let (left_window, right_window) = window_f(max_prob_fold_change);

            let mut diff_exp_distribution = ProbDistribution1d::new();

            let mut calc_prob = |f : f64| {
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


                let density = |i : usize, x :f64 | {
                    prob_dist_i_k1.get(x).ln_add_exp(prob_dist_i_k2.get((f * (x + c)-c)))
                };
    
                let mut prob = LogProb::ln_simpsons_integrate_exp(
                    density,
                    0.,
                    15.,
                    11,
                );

                let prob = LogProb::ln_one();
                diff_exp_distribution.insert(f, prob);

                prob
            };

            // Step 2: For each fold change in determined range, calculate formula 9 and put in ProbabilityDistribution
            for f in left_window {
                calc_prob(f);
            }

            for f in right_window {
                calc_prob(f);
            }

            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;

            Ok(())
        })?;


    Ok(())
}
