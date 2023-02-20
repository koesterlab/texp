//! This implements formula 9 of the document and calculates the fold change / differential expression.
use std::collections::VecDeque;
use std::mem;
use std::path::Path;
use std::fs::File;

use anyhow::Result;
use bio::stats::LogProb;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use csv;
use itertools_num::linspace;
use itertools::iproduct;
use crate::common::{Outdir, Pair, QueryPoints};
use crate::preprocess::Preprocessing;
use crate::prob_distribution_1d::ProbDistribution1d;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::sample_expression;
use std::collections::HashMap;

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

    let query_points = QueryPoints::new(c)?;
    let all_mu_ik_points = query_points.all_mu_ik();
    let start_points_theta_i = query_points.thetas();
    let possible_f = query_points.possible_f();
    let start_points_mu_ik = query_points.start_points_mu_ik();

    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().skip(190400).collect();

    let file = File::open("/vol/nano/bayesian-diff-exp-analysis/texp-evaluation/estimated_dispersion.csv")?;
    let mut rdr = csv::Reader::from_reader(file);
    let mut thetas = Vec::<f64>::new();
    for result in rdr.records() {
        let record = result?;
        let dispersion: f64 = record[1].parse().unwrap();
        thetas.push(dispersion);
    }

    let subsampled_ids = vec!["ERCC-00130","ERCC-00004", "ERCC-00136", "ERCC-00096", "ERCC-00171", "ERCC-00009",
    "ERCC-00074", "ERCC-00113", "ERCC-00145", "ERCC-00002", "ERCC-00046", "ERCC-00003"];

    // feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            // if subsampled_ids.contains(&feature_id.as_str()) {
            // if feature_id.as_str() == "ERCC-00130" { 

            // println!("\n--------------feature {:?} {:?}", i, feature_id);

            // let theta = thetas[i-190432];
            // println!("{:?}", theta);

            let prob_dist_i_k1: ProbDistribution2d = in_dir1.deserialize_value(feature_id)?;
            let prob_dist_i_k2: ProbDistribution2d = in_dir2.deserialize_value(feature_id)?;
            if prob_dist_i_k1.is_na() || prob_dist_i_k2.is_na() {
                println!("skipped {:?}", feature_id);
                return  Ok(());
            }
            // let mut mu1 = 320.;
            // let mut mu2 = 80.;
            // if in_dir1.to_str().unwrap().contains(&"exp_1".to_string()) {
            //     mu1 = 80.;
            //     mu2 = 320.;
            // }

            // let mut start_points_mu_ik = vec![0.];
            // let mut cur_maximum_likelihood_mean = 100. / 20.;
            // if (cur_maximum_likelihood_mean > 0.) {
            //     while cur_maximum_likelihood_mean < 10000. {
            //         start_points_mu_ik.push(cur_maximum_likelihood_mean);
            //         cur_maximum_likelihood_mean = cur_maximum_likelihood_mean * 2.;
            //     }
            // }
            // start_points_mu_ik.push(10000.);
            // println!("start_points_mu_ik {:?}", start_points_mu_ik);
            // let start_points_mu_ik2 = start_points_mu_ik.clone();

            // let mut start_points_theta_i = Vec::<f64>::new();
            // start_points_theta_i.extend(prior.left_window());
            // start_points_theta_i.push(theta);
            // start_points_theta_i.extend(prior.right_window());
            // start_points_theta_i.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // let start_points_theta_i2 = start_points_theta_i.clone();
            // println!("start_points_theta_i {:?}", start_points_theta_i);

            // let normalverteilung = |x: f64, mu: f64| {
            //     let standardabweichung = 10.;
            //     let c = (2.0 * std::f64::consts::PI * standardabweichung * standardabweichung).sqrt();
            //     let c2 = -2.0 * standardabweichung * standardabweichung;
            //     return c * ((x - mu).powi(2) / c2).exp();
            // };

            // let test_calc1 = |mu_ik : f64, theta_i: f64| {
            //     return LogProb::from(normalverteilung(mu_ik, mu1).ln());
            // };

            // let test_calc2 = |mu_ik : f64, theta_i: f64| {
            //     return LogProb::from(normalverteilung(mu_ik, mu2).ln());
            // };


            // let mut prob_dist_i_k1 = ProbDistribution2d::new();
            // prob_dist_i_k1.insert_grid(start_points_mu_ik, start_points_theta_i, test_calc1);
            // println!("prob_dist_i_k1  get 80 theta {:?}", prob_dist_i_k1.get(&[80. ,theta]) );
            // println!("prob_dist_i_k1  get 320 theta{:?}", prob_dist_i_k1.get(&[320. ,theta]) );
            // let mut prob_dist_i_k2 = ProbDistribution2d::new();
            // prob_dist_i_k2.insert_grid(start_points_mu_ik2, start_points_theta_i2, test_calc2);
            // println!("prob_dist_i_k2 max prob {:?}", prob_dist_i_k2.get_max_prob() );
            // println!("prob_dist_i_k2  get 80 theta {:?}", prob_dist_i_k2.get(&[80. ,theta]) );
            // println!("prob_dist_i_k2  get 320 theta{:?}", prob_dist_i_k2.get(&[320. ,theta]) );
            

            let mut prob_d_i_f = ProbDistribution1d::new();
            let mut diff_exp_distribution = ProbDistribution1d::new();

            // let mut start_points_mu_ik: Vec<f64> = linspace(0.1, 1., 10).collect();
            // start_points_mu_ik.extend( linspace(1.5, 100., 198));
            // start_points_mu_ik.extend( linspace(101., 500., 400));
            // // start_points_mu_ik.extend( linspace(600., 3000., 25));
            // // start_points_mu_ik.extend( linspace(4000., 10000., 7));
            // let start_points_mu_ik2 = start_points_mu_ik.clone();
            // let start_points_mu_ik3 = start_points_mu_ik.clone();

            // let mut start_points_theta_i: Vec<f64> = linspace(0.1, 1., 10).collect();
            // start_points_theta_i.extend( linspace(1.5, 10., 18).step_by(2));
            // start_points_theta_i.extend( linspace(11., 165., 155).step_by(10));
            // // println!("start_points_theta_i {:?}", start_points_theta_i);

            // let mut possible_f = Vec::<f64>::new(); 
            // let mut possible_f2 = HashMap::<OrderedFloat<f64>, Vec<f64>>::new();
            // for (mu_1, mu_2) in iproduct!(start_points_mu_ik, start_points_mu_ik2) {
            //     let f = (mu_2+c)/(mu_1+c);
            //     if (f * (mu_1 + c) - c) > 0.{ // && !possible_f.contains(&f) && f < 20. && f> 1./20.
            //         if let Some(list) = possible_f2.get_mut(&OrderedFloat(f)) {
            //             list.push(mu_1);
            //         } else {
            //             possible_f2.insert(OrderedFloat(f), vec![mu_1]);
            //         }
            //         // possible_f.push(f);
            //     }
            // }
            // for (f, list) in possible_f2.clone() {
            //     if list.len() < 10{
            //         possible_f2.remove(&f);
            //     }
            //     possible_f.push(f64::from(f));
            // }
            // possible_f.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic

            // // println!("fs {:?}", possible_f);
            // println!(" len fs {:?}", possible_f2.len());

            // let calc_prob = |f: f64, list_mu| -> LogProb {
            for f in possible_f{
                // let f =f64::from(f);
                let calc_prob_fixed_theta = |theta|{
                    let density_x= |_, x: f64| {
                        // if f >= 0.95 && f <=0.99 { 
                            // println!("f {:?}, x {:?}, f * (x + c) s- c {:?}", f, x, f * (x + c) - c);
                        // }
                        prob_dist_i_k1.get(&[(f * (x + c) - c), theta]) + prob_dist_i_k2.get(&[x, theta])
                    };
                    let prob_x = LogProb::ln_trapezoidal_integrate_grid_exp(
                        density_x,
                        &start_points_mu_ik
                    );
                    prob_x
                };
                
                let density_theta = |_, theta:f64|{
                    calc_prob_fixed_theta(theta) + prior.prob(theta)
                };
                let prob_theta = density_theta(0., 0.2);

                // let prob_theta = LogProb::ln_trapezoidal_integrate_grid_exp(
                //     density_theta,
                //     &start_points_theta_i
                // );
                // prob_theta
            // };
                
                    // let value = calc_prob(f64::from(f), list_mu);
                    // if f >= 0.95 && f <=0.99 { 
                    //     println!("prob_d_i_f f {:?} {:?}", f, value.exp());
                    // }
                    prob_d_i_f.insert(*f, prob_theta);
                }

            let density = |_, f| prob_d_i_f.get(f);

            let prob_f = LogProb::ln_trapezoidal_integrate_grid_exp(
                density,
                &possible_f
            );
            let calc_prob_f = |f| {                
                let prob = prob_d_i_f.get(f) - prob_f;
                prob
            };

            for f in possible_f.clone(){
                let value = calc_prob_f(f);
                // println!("diff_exp_distribution f {:?} {:?}", f, value);
                diff_exp_distribution.insert(f, value );
            }


        //     possible_f.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
        //     let min_f = possible_f[0];
        //     let max_f = possible_f[possible_f.len()-1];
            
        //     possible_f = possible_f.iter().step_by(2500).copied().collect();
        //     println!("len of possible_f {:?}, min {:?}, max {:?}", possible_f.len(),min_f , max_f);
        //     // println!("possible_f {:?}", possible_f);

        //     let calc_prob = |f: f64| {
        //             let calc_prob_fixed_theta = |theta|{
        //                 let density_x= |_, x: f64| {
        //                     // if f >= 0.95 && f <=0.99 { 
        //                         // println!("f {:?}, x {:?}, f * (x + c) s- c {:?}", f, x, f * (x + c) - c);
        //                     // }
        //                     prob_dist_i_k1.get(&[(f * (x + c) - c), theta]) + prob_dist_i_k2.get(&[x, theta])
        //                 };
        //                 // let mut points1 = prob_dist_i_k1.points.keys().map(|value| -c + (c + value.raw()) / f).collect::<Vec<_>>();
        //                 // let mut points2 = prob_dist_i_k2.points.keys().map(|value| value.raw()).collect::<Vec<_>>();
        
        //                 // points1.append(&mut points2);
        //                 // points1.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
        //                 // points1.dedup();              
        
        //                 // let probs_x = points1
        //                 //                 .windows(2)
        //                 //                 .map(|x| LogProb::ln_simpsons_integrate_exp(density_x, x[0], x[1], 5))
        //                 //                 .collect::<Vec<_>>();
        //                 // let prob_x= LogProb::ln_sum_exp(&probs_x);
        //                 // let x_range = prob_dist_i_k1.get_range_per_theta(theta);
        //                 // println!("theta {:?}, x_range {:?}", theta, x_range);
        //                 // let prob_x = LogProb::ln_simpsons_integrate_exp(density_x, 
        //                 //     x_range[0], 
        //                 //     x_range[1], 
        //                 //     301);
                        
        //                 let prob_x = LogProb::ln_trapezoidal_integrate_grid_exp(
        //                     density_x,
        //                     &start_points_mu_ik3
        //                 );
        //                 // if f >= 0.95 && f <=0.99 { 
        //                 //     println!("prob_x {:?}", prob_x.exp());
        //                 // }

        //                 prob_x
        //         };
                
        //         let density_theta = |_, theta:f64|{
        //             calc_prob_fixed_theta(theta) + prior.prob(theta)
        //         };

        //         // let prob_theta = LogProb::ln_simpsons_integrate_exp(
        //         //         density_theta,
        //         //         prior.min_value(),
        //         //         10.,
        //         //         21,
        //         //     );
        //         let prob_theta = density_theta(0., 0.4);

        //         // let prob_theta = LogProb::ln_trapezoidal_integrate_grid_exp(
        //         //     density_theta,
        //         //     &start_points_theta_i
        //         // );


        //         // println!("f {:?}, prob_theta {:?}", f, prob_theta.exp());
        //         // let prob1 = LogProb::ln_simpsons_integrate_exp(
        //         //         density_theta,
        //         //         prior.min_value(),
        //         //         20.,
        //         //         151,
        //         // );
        //         // let prob2 = LogProb::ln_simpsons_integrate_exp(
        //         //     density_theta,
        //         //     20.,
        //         //     prior.max_value(),
        //         //     101,
        //         // );
        //         // let prob = prob1.ln_add_exp(prob2); 
        //         prob_theta
        //     };

        //     for f in possible_f.clone(){
        //             let value = calc_prob(f as f64);
        //             if f >= 0.95 && f <=0.99 { 
        //                 println!("prob_d_i_f f {:?} {:?}", f, value.exp());
        //             }
        //             prob_d_i_f.insert(f as f64, value);
        //     }

        //     let density = |_, f| prob_d_i_f.get(f);

        //     let prob_f = LogProb::ln_trapezoidal_integrate_grid_exp(
        //         density,
        //         &possible_f
        //     );
        //     // let prob_f = LogProb::ln_simpsons_integrate_exp(
        //     //     density,
        //     //     min_f,
        //     //     max_f,
        //     //     51,
        //     // );

        //     // println!("prob_f {:?}", prob_f);

        //     let calc_prob_f = |f| {                
        //         let prob = prob_d_i_f.get(f) - prob_f;
        //         prob
        //     };

        //     for f in possible_f.clone(){
        //         let value = calc_prob_f(f as f64);
        //         // println!("diff_exp_distribution f {:?} {:?}", f, value);
        //         diff_exp_distribution.insert(f as f64, value );
        // }

            // for x in 0..10{
            //     let f = x as f64 / 10. ;
            //     println!("f {:?}", f);
            //     diff_exp_distribution.insert(f as f64, calc_prob_f(f as f64));
            //     diff_exp_distribution.insert(f as f64 + 0.05, calc_prob_f(f as f64 +0.05));
                
            // }
            // for f in 1..11{
            //     println!("f {:?}", f);
            //     diff_exp_distribution.insert(f as f64, calc_prob_f(f as f64));

            //     diff_exp_distribution.insert(f as f64 + 0.25, calc_prob_f(f as f64 +0.25));
            //     diff_exp_distribution.insert(f as f64 + 0.5, calc_prob_f(f as f64 +0.5));
            //     diff_exp_distribution.insert(f as f64 + 0.75, calc_prob_f(f as f64 +0.75));
            // }


            // let mut start_points = vec![0.];
            // let mut cur_prob = LogProb::ln_zero();
            // diff_exp_distribution.insert(0., cur_prob);
            // let mut cur_max_prob_fold_change = 1. / max_prob_fold_change / 16.;
            // if cur_max_prob_fold_change > 0. {
            //     while cur_max_prob_fold_change < 32. {
            //         start_points.push(cur_max_prob_fold_change);
            //         cur_prob = calc_prob(cur_max_prob_fold_change);
            //         diff_exp_distribution.insert(cur_max_prob_fold_change, cur_prob);
            //         // println!("insert mu {:?}, prob {:?}", cur_maximum_likelihood_mean, f64::from(cur_prob.exp()));
            //         cur_max_prob_fold_change = cur_max_prob_fold_change * (2.0_f64).sqrt();
            //     }
            // }
            // if start_points.len() == 1 {
            //     start_points.push(5.);
            //     cur_prob = calc_prob(5.);
            //     diff_exp_distribution.insert(5., cur_prob);
            //     start_points.push(16.);
            //     cur_prob = calc_prob(16.);
            //     diff_exp_distribution.insert(16., cur_prob);
            //     // println!("insert mu {:?}, prob {:?}", 5000., f64::from(cur_prob.exp()));
            // }
            // start_points.push(32.);
            // cur_prob = calc_prob(32.);
            // diff_exp_distribution.insert(32., cur_prob);
            // println!("start points fs {:?}", start_points);


            // let mut queue = VecDeque::<Pair>::new();
            // start_points.windows(2).for_each(|w| {
            //     queue.push_back(Pair {
            //         left: w[0],
            //         right: w[1],
            //     })
            // });

            // while queue.len() > 0 {
            //     // println!("Len of queue {:?}", queue.len());
            //     // println!("#values inserted {:?}", diff_exp_distribution.len());
            //     let pair = queue.pop_front().unwrap();
            //     let left = pair.left;
            //     let right = pair.right;
            //     if left == 0. && right == 0. {
            //         continue;
            //     }
            //     let mut middle = 0.;
            //     if left.is_finite() && right.is_finite() {
            //         middle = left / 2. + right / 2.;
            //     } else if left.is_finite() {
            //         middle = 10. * left;
            //     }
            //     if middle.is_infinite() {
            //         continue;
            //     }
            //     let estimated_value = diff_exp_distribution.get(middle);
            //     let calculated_value = calc_prob(middle);
            //     // println!("middle {:?}, est: {:?}, calc: {:?}", middle, estimated_value.exp(), calculated_value.exp());
            //     diff_exp_distribution.insert(middle, calculated_value);
            //     if (estimated_value.exp() - calculated_value.exp()).abs() > 0.1 {
            //         queue.push_back(Pair {
            //             left: left,
            //             right: middle,
            //         });
            //         queue.push_back(Pair {
            //             left: middle,
            //             right: right,
            //         });
            //     }
            // }
            // in_dir1.serialize_value(feature_id, prob_dist_i_k1)?;
            // in_dir2.serialize_value(feature_id, prob_dist_i_k2)?;

            // Step 3: Write output
            out_dir.serialize_value(feature_id, diff_exp_distribution)?;
        // }
            Ok(())
        })?;

    Ok(())
}
