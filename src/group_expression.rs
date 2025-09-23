//! This implements formula 5,6,7 of the document.
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use csv;
use rayon::prelude::*;
use duckdb::{Connection, params};
use std::sync::{Arc, Mutex};

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
    let sample_ids = preprocessing
        .scale_factors()
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    // let prior = preprocessing.prior()?;
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let db_path = out_dir_path.to_str().unwrap(); //format!("{}.duckdb", out_dir_path.to_str().unwrap());
    let conn = Connection::open(db_path)?;
    ProbDistribution2d::init_schema(&conn)?; // ensure schema exists
    // Wrap in Arc<Mutex<Connection>> for parallel use
    let conn = Arc::new(Mutex::new(conn));

    println!("Before feature_ids par_iter");
    feature_ids
        .par_iter().take(5)
        .try_for_each(|(i, feature_id)| -> Result<()> {
            println!("--------------feature {:?} {:?}", i, feature_id);

            let sample_expression_likelihoods: Vec<ProbDistribution2d> = sample_expression_paths
                .iter()
                .map(|path| {
                    // Each ProbDistribution2d owns its own connection
                    ProbDistribution2d::with_readonly_connection(path.to_str().unwrap(), &feature_id)
                })
                .collect::<duckdb::Result<_>>()?;
            println!("feature {:?} After reading likelihoods", feature_id);
            // let maximum_likelihood_mean = maximum_likelihood_means.iter().sum::<f64>()
            //     / maximum_likelihood_means.len() as f64;

            let prob_dist = ProbDistribution2d::with_connection(conn.clone(), feature_id).unwrap();
            println!("feature {:?} After creating prob_dist", feature_id);

            let calc_prob = |mu_ik: f64, theta_i: f64| {
                // println!("mu_ik {:?}", mu_ik);
                if mu_ik == 0. {
                    return LogProb::ln_zero();
                }
                let prob = sample_expression_likelihoods
                    .iter()
                    .map(|sample_expression_likelihood| {
                        sample_expression_likelihood.get(mu_ik, theta_i) //.unwrap_or(LogProb::ln_zero())
                    })
                    .sum::<LogProb>(); //Formula 5
                //                        // +LogProb(*prior.prob(theta_i));
                //                        // prob = LogProb(f64::from(prob) * 8.);

                // Result of formula 7.
                // let prob= LogProb::ln_simpsons_integrate_exp(
                //     density,
                //     prior.min_value(),
                //     prior.max_value(),
                //     451,
                // );
                // let prob = density(0., prior.mean());
                // println!("mu_ik {:?}, prob {:?}", mu_ik, prob);
                // if theta_i == 0.01 {
                //     wtr.serialize((mu_ik, prob.exp())).unwrap();
                // }
                prob
            };

            let query_points = query_points::calc_query_points(
                c,
                preprocessing.mean_disp_estimates().clone(),
                sample_ids.clone(),
                preprocessing.feature_ids().clone(),
                *i,
            );
            let start_points_mu_ik = query_points.all_mu_ik();
            let start_points_theta_i = query_points.thetas();

            // Compute grid in memory
            println!("feature {:?} before compute_grid mu_ik_points.len() {:?}, start_points_theta_i.len() {:?}", feature_id, start_points_mu_ik.len(), start_points_theta_i.len());
            let probs = prob_dist.compute_grid(&start_points_mu_ik, &start_points_theta_i, calc_prob);
            println!("feature {:?} after compute_grid", feature_id);
            // Write results to DuckDB (mutex ensures serialized access)
            prob_dist.write_output(&probs).unwrap();
            println!("feature {:?} after write_output ", feature_id);

            // let norm_factor = prob_dist.normalize(); // remove factor c_ik
            // out_dir.serialize_value(feature_id, prob_dist)?;
            // }
            Ok(())
        })?;

    Ok(())
}
