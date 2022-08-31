//! This implements formula 3+4 of the document.
use std::fs;
use std::mem;
use std::path::Path;
// use std::collections::VecDeque;
use std::process::exit;


use anyhow::Result;
use bio::stats::LogProb;
use getset::Getters;
use itertools_num::linspace;
// use itertools::iproduct;
use rayon::prelude::*;
use rmp_serde::Deserializer;
use serde::Deserialize as SerdeDeserialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;
// use kdtree::distance::squared_euclidean;

use crate::common::{window_x, Outdir}; // , Square, Point
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::ProbDistribution2d;

pub(crate) fn sample_expression(
    preprocessing: &Path,
    sample_id: &str,
    epsilon: LogProb,
    out_dir: &Path,
) -> Result<()> {
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let prior = preprocessing.prior()?;
    let mean_disp_estimates =
        preprocessing
            .mean_disp_estimates()
            .get(sample_id)
            .ok_or(Error::UnknownSampleId {
                sample_id: sample_id.to_owned(),
            })?;
    let group_means = preprocessing.group_means();

    let s_j = preprocessing
        .scale_factors()
        .get(sample_id)
        .unwrap()
        .clone();

    let mut feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    let out_dir = Outdir::create(out_dir)?;
    feature_ids.truncate(10000);
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            if feature_id.as_str() == "ENST00000533357.5" {

            let d_ij = mean_disp_estimates.means()[*i]; //TODO Do we need group mean mu_ik instead of sample mean
            let mut d_ik = group_means[*i];
            if d_ik < 1e-8 { //TODO SEnsible? Filter at different position better?
                d_ik = 0.;
            }
            // METHOD: If the per-sample dispersion is unknown, fall back to a mean interpolated from the other samples.
            let t_ij = if let Some(t_ij) = mean_disp_estimates.dispersions()[*i] {
                t_ij
            } else if let Some(t_ij) = preprocessing.interpolate_dispersion(*i) {
                t_ij
            } else {
                // TODO log message
                return Ok(());
            };
            println!("\n\n -------------------------- feature {:?}", i);
            println!("d_ij {:?}, d_ik {:?}, t_ij {:?}", d_ij, d_ik, t_ij);

            // if *i == 246_usize {

            // let max_prob = prob_mu_ik_theta_i_x(d_ij, d_ij, d_ij, t_ij, prior.mean(), s_j);
            // let prob_threshold = LogProb(*max_prob - 10.0f64.ln());
            // let (mu_ik_left_window, mu_ik_right_window) = window(d_ik);

            let mut likelihoods = ProbDistribution2d::new();
            // println!("--------------------------3");
            let mut start_points_mu_ik = vec![
                0.,
                // d_ik / 5.,
                // d_ik,
                // d_ik * 10.,
                // d_ik * 100.,
            ];
            if d_ik != 0. {
                let mut cur_d_ik = d_ik / 10.;
                // for (int i = 0; i < 4 && cur_d_ik < 10000; ++i, cur_d_ik *= 10.)
                for i in 1..4  {
                    if cur_d_ik >= 10000. { break; }
                    start_points_mu_ik.push(cur_d_ik);
                    cur_d_ik = cur_d_ik * 10.;
                }
                if start_points_mu_ik.len() == 1 {
                    start_points_mu_ik.push(5000.);
                }
                start_points_mu_ik.push(10000.);
            }
            println!("start_points_mu_ik {:?}",start_points_mu_ik);
            let mut start_points_theta_i = Vec::<f64>::new();
            start_points_theta_i.extend(prior.left_window());
            start_points_theta_i.extend(prior.right_window());
            start_points_theta_i.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // println!("start_points_theta_i {:?}", start_points_theta_i);
            // println!("--------------------------4");
            // for (m, t) in iproduct!(&start_points_mu_ik, &start_points_theta_i) {
            //     let prob = likelihood_mu_ik_theta_i(
            //         d_ij,
            //         *m,  // mu_ik
            //         t_ij,
            //         *t,  // *theta_i
            //         s_j,
            //         epsilon,
            //     );
            //     println!("mu {:?}, theta {:?}, prob {:?}", m ,t ,prob);
            //     likelihoods.insert(*m, *t, prob);
            // }

            let calc_prob = |m, t| {
                likelihood_mu_ik_theta_i(
                    d_ij,
                    m,  // mu_ik
                    t_ij,
                    t,  // theta_i
                    s_j,
                    epsilon,
                )
            };
            likelihoods.insert_grid(start_points_mu_ik, start_points_theta_i, calc_prob);

            // start_points_mu_ik.push(f64::INFINITY);
            // start_points_theta_i.push(f64::INFINITY);

            // let mut queue = VecDeque::<Square>::new();
            // for (m, t) in iproduct!(start_points_mu_ik.windows(2), start_points_theta_i.windows(2)) {
            //     queue.push_back(Square { top_left: Point{x: m[0], y: t[1]},
            //         top_right: Point{x: m[1], y: t[1]},
            //         bot_left: Point{x: m[0], y: t[0]},
            //         bot_right: Point{x: m[1], y: t[0]},
            //     });
            // }

            // while queue.len() > 0 {
            //     println!("\n\n len {:?}", queue.len());
            //     // println!("--------------------------");
            //     let square = queue.pop_front().unwrap();
            //     // println!("square tl {:?}, tr {:?}, bl {:?}, br {:?}", square.top_left, square.top_right, square.bot_left, square.bot_right);
            //     println!("square {:?}", square);
            //     let new_point_candidates = new_points(&square);
            //     let mut calced_values = Vec::<Option<LogProb>>::new();
            //     let mut estimated_values = Vec::<Option<LogProb>>::new();
            //     for p in &new_point_candidates {
            //         if p.x.is_infinite() || p.y.is_infinite(){
            //             calced_values.push(None);
            //             estimated_values.push(None);
            //             continue;
            //         }
            //         let distance_to_nearest = likelihoods.kdtree.nearest(&[p.x, p.y], 1, &squared_euclidean).unwrap()[0].0;
            //         println!("p.x {:?}, p.y {:?}, dist {:?}", p.x, p.y, distance_to_nearest);
            //         if distance_to_nearest / squared_euclidean(&[p.x, p.y], &[square.bot_left.x, square.bot_left.y]) > 1e-6 {
            //             let prob = likelihood_mu_ik_theta_i(
            //                 d_ij,
            //                 p.x,  // mu_ik
            //                 t_ij,
            //                 p.y,  // *theta_i
            //                 s_j,
            //                 epsilon,
            //             );
            //             // println!("prob {:?}, d_ij {:?}, p.x {:?}, t_ij {:?}, p.y {:?}", prob, d_ij, p.x, t_ij, p.y);
            //             calced_values.push(Some(prob));
            //             estimated_values.push(Some(likelihoods.get(&[p.x, p.y])));
            //         } else {
            //             calced_values.push(None);
            //             estimated_values.push(None);
            //         }
            //     }
            //     let mut need_iteration = false;
            //     for (estimated_value, calculated_value) in estimated_values.iter().zip(calced_values.iter()) {
            //         if estimated_value.is_none() {
            //             continue;
            //         }
            //         let prob_dist_max = likelihoods.get_max_prob();
            //         if difference_to_big(estimated_value.unwrap(), calculated_value.unwrap(), prob_dist_max) {
            //             need_iteration = true;
            //             break;
            //         }
            //     }
            //     // // println!("middle {:?},est {:?}, calc {:?}, len {:?}", middle, estimated_value, calculated_value, queue.len());
            //     for (point, calculated_value) in new_point_candidates.iter().zip(calced_values.iter()) {
            //         if calculated_value.is_none() {
            //             continue;
            //         }
            //         likelihoods.insert(point.x, point.y, calculated_value.unwrap());
            //     }
            //     if need_iteration {
            //         queue.push_back(Square { top_left: square.top_left,
            //             top_right: new_point_candidates[4],
            //             bot_left: new_point_candidates[1],
            //             bot_right: new_point_candidates[2],
            //         });
            //         queue.push_back(Square { top_left: new_point_candidates[4],
            //             top_right: square.top_right,
            //             bot_left: new_point_candidates[2],
            //             bot_right: new_point_candidates[3],
            //         });
            //         queue.push_back(Square { top_left: new_point_candidates[1],
            //             top_right: new_point_candidates[2],
            //             bot_left: square.bot_left,
            //             bot_right: new_point_candidates[0],
            //         });
            //         queue.push_back(Square { top_left: new_point_candidates[2],
            //             top_right: new_point_candidates[3],
            //             bot_left: new_point_candidates[0],
            //             bot_right: square.bot_right,
            //         });
            //     }
            // }

            out_dir.serialize_value(feature_id, likelihoods)?;
            }
            Ok(())
            
        })?;

    out_dir.serialize_value(
        "info",
        SampleInfo {
            sample_id: sample_id.to_owned(),
        },
    )?;

    Ok(())
}

// #[allow(unused)]
// pub(crate) fn new_points(
//     square: &Square,
// ) -> Vec<Point> {
//     let mut result = Vec::<Point>::new();

//     let new_point = |p1: Point, p2: Point| {
//         let new_coord = |left: f64, right: f64| {
//             let middle; //= 0.;
//             if right.is_finite() {
//                 middle = left / 2. + right / 2.;
//             } else {
//                 middle = 10. * left;
//             }
//             return middle;
//         };
//         return Point{x:new_coord(p1.x, p2.x), y:new_coord(p1.y, p2.y)};
//     };

//     result.push(new_point(square.bot_left, square.bot_right));
//     result.push(new_point(square.bot_left, square.top_left));
//     result.push(new_point(square.bot_left, square.top_right));
//     result.push(new_point(square.bot_right, square.top_right));
//     result.push(new_point(square.top_left, square.top_right));

//     return result;
// }

#[derive(Debug, Deserialize, Serialize, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct SampleInfo {
    sample_id: String,
}

impl SampleInfo {
    #[allow(unused)]
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(SampleInfo::deserialize(&mut Deserializer::new(
            fs::File::open(path)?,
        ))?)
    }
}

fn prob_mu_ik_theta_i_x(
    x: f64,
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    s_j: f64,
) -> LogProb {
    if d_ij == 0. {
        println!("############ prob_mu_ik_theta_i_x x {:?}, d_ij {:?}, mu_ik {:?}, t_ij {:?}, theta_i {:?}, s_j {:?}", x, d_ij, mu_ik, t_ij, theta_i, s_j);
        println!("neg_binom(d_ij, x, t_ij) {:?}, neg_binom(x, mu_ik * s_j, theta_i) {:?},result {:?}", neg_binom(d_ij, x, t_ij), neg_binom(x, mu_ik * s_j, theta_i),neg_binom(d_ij, x, t_ij) + neg_binom(x, mu_ik * s_j, theta_i));
    }
    neg_binom(d_ij, x, t_ij) + neg_binom(x, mu_ik * s_j, theta_i)
}

/// Inner of equation 3/4 in the document.
fn likelihood_mu_ik_theta_i(
    d_ij: f64,
    mu_ik: f64,
    t_ij: f64,
    theta_i: f64,
    s_j: f64,
    _: LogProb,   // epsilon
) -> LogProb {
    if mu_ik == 0. {
        return LogProb::ln_zero();
    }
    if mu_ik == 3208600050.4032364 {
        println!("############ likelihood_mu_ik_theta_i");
        println!("mu_ik {:?}, theta_i {:?}, d_ij {:?}, t_ij {:?}, s_j {:?}", mu_ik, theta_i, d_ij, t_ij, s_j);
    }
    // TODO determine whether this is the best window given that we also have access to mu_ik here.
    // let (x_left_window, x_right_window) = window_x(d_ij);
    // let window = linspace::<f64>(0., 1001., 1000);
    let mut max_prob = LogProb::ln_zero();
    // let prob = |x| {
    //     let calced_prob = prob_mu_ik_theta_i_x(x, d_ij, mu_ik, t_ij, theta_i, s_j);
    //     if calced_prob > max_prob{
    //         max_prob = calced_prob
    //     }       
    //     calced_prob
    // }; 

    // let max_prob = prob(d_ij);
    // let threshold = LogProb(*max_prob - 10.0f64.ln());
    // let is_informative = |prob: &LogProb| *prob >= threshold; // TODO maybe better relative to the maximum?

    let mut result = LogProb::ln_zero();
    let mut x : u64 = 0;
    while true{
        let calced_prob = prob_mu_ik_theta_i_x(x as f64, d_ij, mu_ik, t_ij, theta_i, s_j);
        if calced_prob > max_prob{
            max_prob = calced_prob;
        }       
        result = result.ln_add_exp(calced_prob);
        if x > 500 && calced_prob - max_prob < LogProb(0.01_f64.ln()){
            break;
        }
        // if x as i64 % 1000 == 0{
        //     println!("x {:?}, calced prob {:?}, max_prob {:?}", x, calced_prob, max_prob);
        //     println!("mu_ik {:?}, theta_i {:?}, d_ij {:?}, t_ij {:?}, s_j {:?}", mu_ik, theta_i, d_ij, t_ij, s_j);
        // } 
        x = x+1;
    }
    // println!("mu_ik {:?}, theta_i {:?}, x {:?}", mu_ik, theta_i, x);
    result
    // LogProb::ln_sum_exp(&probs)
    // LogProb::ln_sum_exp(
    //     &window
    //         .map(&prob)
    //         // .take_while(&is_informative)
    //         .collect::<Vec<_>>(),
    // )
}

fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);

    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    // (p1 - b + p2).exp() / (x + n)
    if x == 0. {
        // println!("x {:?}, mu {:?}, theta {:?}, n {:?}, p {:?}, p1 {:?}, p2 {:?}, b {:?}", x, mu, theta, n, p ,p1, p2, b);
    }
    LogProb((p1 - b + p2) - (x + n).ln()) // TODO is this the correct form for returning LogProb?
}

// def neg_binom(x, mu, theta):
//     """ Own implementation of the negative negative binomial distribution using betaln
//     """
//     n = 1. / theta
//     p = n / (n + mu)
//     p1 = n*np.log(p) if n > 0 else 0
//     p2 = x*np.log(1-p) if x > 0 else 0
//     b = betaln(x + 1, n)
//     if (p1 < p2):
//         return exp(p2 - b + p1) / (x+n)
//     else:
//         return exp(p1 - b + p2) / (x+n)
