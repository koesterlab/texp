//! This implements formula 5,6,7 of the document.
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use rayon::prelude::*;
use csv;

use crate::common::Outdir;
// use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::ProbDistribution2d;
// use crate::sample_expression::SampleInfo;
use crate::query_points;



pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    c: f64,
    out_dir_path: &Path,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let sample_ids = preprocessing.scale_factors().keys().cloned().collect::<Vec<_>>();
    // let prior = preprocessing.prior()?;
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    // let file = File::open("/vol/nano/bayesian-diff-exp-analysis/texp-evaluation/estimated_dispersion.csv")?;
    // let mut rdr = csv::Reader::from_reader(file);
    // let mut thetas = Vec::<f64>::new();
    // for result in rdr.records() {
    //     let record = result?;
    //     let dispersion: f64 = record[1].parse().unwrap();
    //     thetas.push(dispersion);
    // }

    let out_dir = Outdir::create(out_dir_path)?;
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {

            // println!("--------------feature {:?} {:?}", i, feature_id);
            // let maximum_likelihood_means: Vec<f64> = sample_expression_paths
            //     .iter()
            //     .map(|sample_expression_path| {
            //         let sample_info: SampleInfo =
            //             Outdir::open(sample_expression_path)?.deserialize_value("info")?;
            //         Ok(preprocessing
            //             .mean_disp_estimates()
            //             .get(sample_info.sample_id())
            //             .ok_or(Error::UnknownSampleId {
            //                 sample_id: sample_info.sample_id().to_owned(),
            //             })?
            //             .means()[*i])
            //     })
            //     .collect::<Result<Vec<_>>>()?;
            // println!("1");
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
                if sample_expression_likelihoods.iter().all(|x| x.is_na()){
                    out_dir.serialize_value(feature_id, ProbDistribution2d::na())?;
                    return Ok(());
                }

            // let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
            //     / maximum_likelihood_means.len() as f64;

            let mut prob_dist = ProbDistribution2d::new();
            // Extend out_dir_path with feature_id and extension csv
            let mut output = out_dir_path.to_path_buf();
            output.push(feature_id);
            output.set_extension("csv");

            let mut wtr = csv::Writer::from_path(output)?;
            wtr.serialize(("mu_ik", "probability")).unwrap();
            let calc_prob = |mu_ik : f64, theta_i: f64| {
                // println!("mu_ik {:?}", mu_ik);
                if mu_ik == 0. {
                    return LogProb::ln_zero();
                }
                let prob = sample_expression_likelihoods
                    .iter()
                    .map(|sample_expression_likelihood| {
                        sample_expression_likelihood
                        .get(&[mu_ik, theta_i])
                    })
                    .sum::<LogProb>();  //Formula 5
                    // +LogProb(*prior.prob(theta_i));
                // prob = LogProb(f64::from(prob) * 8.);

                // Result of formula 7.
                // let prob= LogProb::ln_simpsons_integrate_exp(
                //     density,
                //     prior.min_value(),
                //     prior.max_value(),
                //     451,
                // );
                // let prob = density(0., prior.mean());
                // println!("mu_ik {:?}, prob {:?}", mu_ik, prob);
                if theta_i == 0.01 {
                    wtr.serialize((mu_ik, prob.exp())).unwrap();
                }
                prob

            };

            let query_points = query_points::calc_query_points(c, preprocessing.mean_disp_estimates().clone(), sample_ids.clone(), preprocessing.feature_ids().clone(), *i);
            let start_points_mu_ik = query_points.all_mu_ik();
            let start_points_theta_i = query_points.thetas();

            prob_dist.insert_grid(start_points_mu_ik.clone(), start_points_theta_i.clone(), calc_prob);

            // let norm_factor = prob_dist.normalize(); // remove factor c_ik
            out_dir.serialize_value(feature_id, prob_dist)?;
            // }
            Ok(())
        })?;

    Ok(())
}
