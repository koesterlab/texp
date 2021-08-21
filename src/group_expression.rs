//! This implements formula 5,6,7 of the document.

use std::path::{Path, PathBuf};
use std::fs;

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;

use crate::common::{ProbDistribution, window};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prior::Prior;
use crate::sample_expression::SampleExpression;

pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    prior: &Prior,
    out_dir: &Path,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let feature_ids = preprocessing.feature_ids();

    let sample_expressions: Vec<SampleExpression> = sample_expression_paths
        .iter()
        .map(|path| SampleExpression::from_path(path))
        .collect::<Result<Vec<_>>>()?;
    
    if out_dir.exists() {
        return Err(Error::ExistingOutputDir { path: out_dir.to_owned() }.into());
    }
    fs::create_dir_all(out_dir)?;

    for (i, feature_id) in feature_ids.iter().enumerate() {
        let maximum_likelihood_means: Vec<f64> = sample_expressions
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
            .collect::<Result<Vec<_>>>()?;

        let maximum_likelihood_mean =
            maximum_likelihood_means.iter().sum::<f64>() / maximum_likelihood_means.len() as f64;

        let (left_window, right_window) = window(maximum_likelihood_mean);

        let mut prob_dist = ProbDistribution::default();

        let mut calc_prob = |mu_ik| {
            let density = |i, theta_i| {
                sample_expressions
                    .iter()
                    .map(|sample_expression: &SampleExpression| {
                        let likelihood_func: &ProbDistribution<(N32, N32)> = sample_expression.likelihoods().get(i).unwrap();
                        likelihood_func.get(&(N32::new(mu_ik as f32), N32::new(theta_i as f32))) + prior.prob(theta_i)
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

            prob_dist.insert(N32::new(mu_ik as f32), prob);

            prob
        };

        for mu_ik in left_window {
            calc_prob(mu_ik);
        }

        for mu_ik in right_window {
            calc_prob(mu_ik);
        }

        let file = fs::File::create(out_dir.join(feature_id).with_extension("mpk"))?;
        prob_dist.serialize(&mut Serializer::new(file))?;
    }

    Ok(())
}
