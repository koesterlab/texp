//! This implements formula 5,6,7 of the document.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rayon::prelude::*;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;

use crate::common::{window, Outdir, ProbDistribution};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prior::Prior;
use crate::sample_expression::SampleInfo;

pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    out_dir: &Path,
) -> Result<()> {

    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let out_dir = Outdir::create(out_dir)?;

    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
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
                    let fullpath = format!("{}{}", sample_expression_path.to_str().unwrap(), feature_id_with_mpk);
                    if (Path::new(&fullpath).exists()) {
                       let likelihoods: ProbDistribution<(N32, N32)> =
                           dir.deserialize_value(feature_id)?;
                       Ok(likelihoods)
                    } else {
                        // TODO Sensible handling of skipped features
                        println!("Skipping {}", fullpath);
                        Ok(ProbDistribution::default())
                    }

                })
                .collect::<Result<Vec<_>>>()?;
            let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
                / maximum_likelihood_means.len() as f64;

            let (left_window, right_window) = window(maximum_likelihood_mean);

            let mut prob_dist = ProbDistribution::default();

            let mut calc_prob = |mu_ik| {
                let density = |i, theta_i| {
                    sample_expression_likelihoods
                        .iter()
                        .map(|sample_expression_likelihood| {
                            sample_expression_likelihood
                                .get(&(N32::new(mu_ik as f32), N32::new(theta_i as f32)))
                                + prior.prob(theta_i)
                        })
                        .sum()  //Formula 5
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
            // TODO ensure that this is a posterior by dividing with integral of all values.

            out_dir.serialize_value(feature_id, prob_dist)?;

            Ok(())
        })?;

    Ok(())
}
