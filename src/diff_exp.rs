//! This implements formula 9 of the document and calculates the fold change / differential expression.

use std::path::Path;

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rayon::prelude::*;

use crate::common::{window, Log2FoldChange, Mean, Outdir, ProbDistribution};
use crate::preprocess::Preprocessing;

pub(crate) fn diff_exp(
    c: f32,
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
            let prob_dist_i_k1: ProbDistribution<Mean> = in_dir1.deserialize_value(feature_id)?;
            let prob_dist_i_k2: ProbDistribution<Mean> = in_dir2.deserialize_value(feature_id)?;

            let mut max_prob_value2 = *prob_dist_i_k2.get_max_prob_value().unwrap();
            let max_prob_value2 = if max_prob_value2 != 0. {
                max_prob_value2
            } else {
                N32::new(1.)
            }; // Avoid division by zero
            let max_prob_fold_change =
                *prob_dist_i_k1.get_max_prob_value().unwrap() / max_prob_value2;

            // Step 1: use window() to determine range around max_prob_fold_change
            let (left_window, right_window) = window(f64::from(max_prob_fold_change));

            let mut diff_exp_distribution = ProbDistribution::<Log2FoldChange>::default();

            let mut calc_prob = |f| {
                 // for x, iterate over values in prob_dist_i_k1, the value for looking up in prob_dist_i_k2 is f * (x + c) - c
                let prob = LogProb::ln_trapezoidal_integrate_grid_exp(
                    |i, value| {
                        prob_dist_i_k1.get(&Mean::new(value))
                            + prob_dist_i_k2.get(&Mean::new(
                                (N32::new(f as f32) * (value + N32::new(c)))
                                    - N32::new(c),
                            ))
                    },
                    &prob_dist_i_k1
                        .points
                        .keys()
                        .map(|value| **value)
                        .collect::<Vec<_>>(),
                );

                diff_exp_distribution.insert(Log2FoldChange::new(N32::new(f as f32)), prob);

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
