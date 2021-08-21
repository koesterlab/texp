//! This implements formula 5,6,7 of the document.

use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;

use crate::common::window;
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prior::Prior;
use crate::sample_expression::SampleExpression;

pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    prior: &Prior,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let feature_ids = preprocessing.feature_ids();

    let sample_expressions: Vec<_> = sample_expression_paths
        .iter()
        .map(|path| SampleExpression::from_path(path))
        .collect()?;

    for i in 0..feature_ids.len() {
        let maximum_likelihood_means: Vec<_> = sample_expressions
            .iter()
            .map(|sample_expression: &SampleExpression| {
                let sample_id = sample_expression.sample_id();
                Ok(preprocessing
                    .mean_disp_estimates()
                    .get(sample_expression.sample_id())
                    .ok_or(Error::UnknownSampleId {
                        sample_id: sample_id.to_owned(),
                    })?
                    .means()[i])
            })
            .collect()?;

        let maximum_likelihood_mean =
            maximum_likelihood_means.iter().sum() / maximum_likelihood_means.len();

        let (left_window, right_window) = window(maximum_likelihood_mean);

        for mu_ik in left_window {
            let density = |i, theta_i| {
                sample_expressions
                    .iter()
                    .map(|sample_expression| {
                        sample_expression.likelihoods()[i].get(mu_ik).get(theta_i)
                    })
                    .sum()
            };

            // Result of formula 7.
            let prob = LogProb::ln_simpsons_integrate_exp(
                density,
                prior.min_value(),
                prior.max_value(),
                11,
            );
        }

        // TODO right window
    }

    Ok(())
}
