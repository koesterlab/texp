//! This implements formula 9 of the document and calculates the fold change / differential expression.
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;
use anyhow::Result;
use bio::stats::LogProb;
use duckdb::Connection;
use noisy_float::types::N64;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

pub(crate) fn diff_exp(
    c: f64,
    preprocessing: &Path,
    group_path1: &Path,
    group_path2: &Path,
    out_dir: &Path,
) -> Result<()> {
    // let out_dir = Outdir::create(out_dir)?;
    let db_path = out_dir.to_str().unwrap(); //format!("{}.duckdb", out_dir_path.to_str().unwrap());
    let conn = Connection::open(db_path)?;
    ProbDistribution1d::init_schema(&conn)?; // ensure schema exists
    // Wrap in Arc<Mutex<Connection>> for parallel use
    let conn = Arc::new(Mutex::new(conn));

    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let _sample_ids = preprocessing
        .scale_factors()
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    let prior = preprocessing.prior()?;
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let query_points_per_feature = query_points::QueryPointsPerFeature::new(&preprocessing, c);

    println!("Before feature_ids par_iter");
    feature_ids
        // .par_iter()
        // .try_for_each(|(i, feature_id)| -> Result<()> {
        .par_chunks(10)
        .try_for_each(|chunk| -> Result<()> {
            for (i, feature_id) in chunk {
                let prob_dist_i_k1_db = ProbDistribution2d::with_readonly_connection(
                    group_path1.to_str().unwrap(),
                    feature_id,
                )?;
                let prob_dist_i_k2_db = ProbDistribution2d::with_readonly_connection(
                    group_path2.to_str().unwrap(),
                    feature_id,
                )?;
                let prob_dist_i_k1 = prob_dist_i_k1_db.load_lookup_table()?;
                let prob_dist_i_k2 = prob_dist_i_k2_db.load_lookup_table()?;

                // if prob_dist_i_k1_db.is_na() || prob_dist_i_k2_db.is_na() {
                //     println!("skipped {:?}", feature_id);
                //     return Ok(());
                // }

                let query_points = query_points_per_feature.get(*i);
                let possible_f = query_points.possible_f();
                let start_points_mu_ik = query_points.start_points_mu_ik();
                let start_points_theta_i = query_points.thetas();

                let mut prob_d_i_f = BTreeMap::<N64, LogProb>::new();
                let diff_exp_distribution =
                    ProbDistribution1d::with_connection(conn.clone(), feature_id).unwrap();

                // let calc_prob = |f: f64, list_mu| -> LogProb {
                for f in possible_f.clone().iter() {
                    // let f =f64::from(f);
                    let calc_prob_fixed_theta = |theta| {
                        let density_x = |_, x: f64| {
                            let mut fx = f * (x + c) - c;
                            if fx < 0.1 {
                                // round fx to 3 decimals
                                fx = (fx * 1000.).round() / 1000.;
                            } else if fx < 100. {
                                // round fx to 2 decimals
                                fx = (fx * 100.).round() / 100.;
                            } else {
                                // round fx to 1 decimal
                                fx = (fx * 10.).round() / 10.;
                            }

                            let p1 = prob_dist_i_k1
                                .get(&(OrderedFloat(fx), OrderedFloat(theta)))
                                .cloned()
                                .unwrap_or(LogProb::ln_zero());
                            let p2 = prob_dist_i_k2
                                .get(&(OrderedFloat(x), OrderedFloat(theta)))
                                .cloned()
                                .unwrap_or(LogProb::ln_zero());

                            p1 + p2
                        };

                        LogProb::ln_trapezoidal_integrate_grid_exp(density_x, start_points_mu_ik)
                    };

                    let density_theta =
                        |_, theta: f64| calc_prob_fixed_theta(theta) + prior.prob(theta);
                    // let prob_theta = density_theta(0., 0.01);

                    let prob_theta = LogProb::ln_trapezoidal_integrate_grid_exp(
                        density_theta,
                        start_points_theta_i,
                    );
                    // prob_theta
                    // };

                    // let value = calc_prob(f64::from(f), list_mu);
                    let f = N64::new(*f);
                    prob_d_i_f.insert(f, prob_theta);
                }

                let density = |_, f| *prob_d_i_f.get(&N64::new(f)).unwrap();

                let prob_f = LogProb::ln_trapezoidal_integrate_grid_exp(density, possible_f);
                let calc_prob_f = |f| {
                    let noisy_f = N64::new(f);
                    let prob = prob_d_i_f.get(&noisy_f).unwrap() - prob_f;
                    prob
                };

                for f in possible_f.clone() {
                    let value = calc_prob_f(f);
                    // println!("feature_id {:?} diff_exp_distribution f {:?} {:?}", feature_id, f, value);
                    diff_exp_distribution.insert(f, value)?;
                }
            }
            Ok(())
        })?;

    Ok(())
}
