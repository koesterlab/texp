//! This implements formula 5,6,7 of the document.

// use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rayon::prelude::*;

use crate::common::{window, Mean, MeanDispersionPair, Outdir};
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
                        // println!("{:?}", feature_id);
                        Ok(ProbDistribution2d::na())
                    }
                })
                .collect::<Result<Vec<_>>>()?;
            let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
                / maximum_likelihood_means.len() as f64;

            let (left_window, right_window) = window(maximum_likelihood_mean);

            let mut prob_dist = ProbDistribution1d::new();

            let mut calc_prob = |mu_ik| {
                let density = |i, theta_i| {
                    let d = sample_expression_likelihoods
                        .iter()
                        .map(|sample_expression_likelihood| {
                            sample_expression_likelihood
                                .get(&[mu_ik, theta_i])
                        })
                        .sum::<LogProb>() + //Formula 5
                        LogProb(*prior.prob(theta_i) * 2.0); // square of Pr(theta_i), formula 8
                        // println!("i {:?}, mu_ik {:?}, theta_i {:?}, summand_theta_i {:?}, density {:?}",i, mu_ik, theta_i, LogProb(*prior.prob(theta_i) * 2.0),  d);
                    d
                };

                // Result of formula 7.
                let prob = LogProb::ln_simpsons_integrate_exp(
                    density,
                    prior.min_value(),
                    prior.max_value(),
                    11,
                );
                prob_dist.insert(mu_ik, prob);

                prob
            };

            // println!("\n\n\nfeature id: {:?} ", feature_id);
            for mu_ik in left_window {
                calc_prob(mu_ik);
            }

            for mu_ik in right_window {
                calc_prob(mu_ik);
            }

            let norm_factor = prob_dist.normalize(); // remove factor c_ik
            if **feature_id == String::from("ENST00000671775.2") {
                println!("norm faktor {:?}", norm_factor);
            }
            out_dir.serialize_value(feature_id, prob_dist)?;
            Ok(())
        })?;

    Ok(())
}
