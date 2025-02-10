//! This implements formula 9 of the document and calculates the fold change / differential expression.
use anyhow::Result;
use bio::stats::LogProb;
use rayon::prelude::*;
use std::path::Path;

use crate::common::Outdir;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;

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
    let sample_ids = preprocessing
        .scale_factors()
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    let prior = preprocessing.prior()?;
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    // let file = File::open("/vol/nano/bayesian-diff-exp-analysis/texp-evaluation/estimated_dispersion.csv")?;
    // let mut rdr = csv::Reader::from_reader(file);
    // let mut thetas = Vec::<f64>::new();
    // for result in rdr.records() {
    //     let record = result?;
    //     let dispersion: f64 = record[1].parse().unwrap();
    //     thetas.push(dispersion);
    // }

    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            // println!("\n--------------feature {:?} {:?}", i, feature_id);

            let prob_dist_i_k1: ProbDistribution2d = in_dir1.deserialize_value(feature_id)?;
            let prob_dist_i_k2: ProbDistribution2d = in_dir2.deserialize_value(feature_id)?;

            if prob_dist_i_k1.is_na() || prob_dist_i_k2.is_na() {
                println!("skipped {:?}", feature_id);
                return Ok(());
            }
            // println!("max_prob keys 1 {:?}", prob_dist_i_k1.get_max_prob_keys());
            // println!("max_prob_1 {:?}", prob_dist_i_k1.get_max_prob().exp());
            // println!("max_prob keys 2 {:?}", prob_dist_i_k2.get_max_prob_keys());
            // println!("max_prob_2 {:?}", prob_dist_i_k2.get_max_prob().exp());
            let query_points = query_points::calc_query_points(
                c,
                preprocessing.mean_disp_estimates().clone(),
                sample_ids.clone(),
                preprocessing.feature_ids().clone(),
                *i,
            );
            let possible_f = query_points.possible_f();
            let start_points_mu_ik = query_points.start_points_mu_ik();
            let start_points_theta_i = query_points.thetas();

            let mut prob_d_i_f = ProbDistribution1d::new();
            let mut diff_exp_distribution = ProbDistribution1d::new();

            // let calc_prob = |f: f64, list_mu| -> LogProb {
            for f in possible_f.clone() {
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

                        let p1 = prob_dist_i_k1.get(&[fx, theta]);
                        let p2 = prob_dist_i_k2.get(&[x, theta]);
                        p1 + p2
                    };
                    let prob_x =
                        LogProb::ln_trapezoidal_integrate_grid_exp(density_x, &start_points_mu_ik);
                    prob_x
                };

                let density_theta =
                    |_, theta: f64| calc_prob_fixed_theta(theta) + prior.prob(theta);
                // let prob_theta = density_theta(0., 0.01);

                let prob_theta = LogProb::ln_trapezoidal_integrate_grid_exp(
                    density_theta,
                    &start_points_theta_i,
                );
                // prob_theta
                // };

                // let value = calc_prob(f64::from(f), list_mu);
                prob_d_i_f.insert(f, prob_theta);
            }

            let density = |_, f| prob_d_i_f.get(f);

            let prob_f = LogProb::ln_trapezoidal_integrate_grid_exp(density, &possible_f);
            let calc_prob_f = |f| {
                let prob = prob_d_i_f.get(f) - prob_f;
                prob
            };

            for f in possible_f.clone() {
                let value = calc_prob_f(f);
                // println!("feature_id {:?} diff_exp_distribution f {:?} {:?}", feature_id, f, value);
                diff_exp_distribution.insert(f, value);
            }

            // println!("prob_d_i_f get_max_prob_position {:?}", prob_d_i_f.get_max_prob_position());
            // println!("diff exp get_max_prob_position {:?}", diff_exp_distribution.get_max_prob_position());

            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;
            // }
            Ok(())
        })?;

    Ok(())
}
