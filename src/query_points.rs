use std::collections::HashMap;

use anyhow::Result;
use getset::Getters;
use itertools_num::linspace;
use itertools::iproduct;
use serde_derive::{Deserialize, Serialize};
use ndarray::{Array1, Dim};

use crate::preprocess::Estimates;


#[derive(Serialize, Deserialize, Debug, Getters)]
pub(crate) struct QueryPoints {
    #[get = "pub(crate)"]
    start_points_mu_ik: Vec<f64> ,
    #[get = "pub(crate)"]
    all_mu_ik: Vec<f64> ,
    #[get = "pub(crate)"]
    thetas: Vec<f64>,
    #[get = "pub(crate)"]
    possible_f: Vec<f64>,
}

impl QueryPoints {
    // get query points for one feature
    pub(crate) fn new(c: f64, mean_mu:f64, min_max: (f64,f64)) -> Result<Self> {
        let mut min = min_max.0;
        let max = min_max.1;
        let mut start_points_mu_ik =  vec![0.];
        if mean_mu < 0.1{
            start_points_mu_ik.extend( linspace(0.001, 0.2, 50).step_by(1));
        } else if mean_mu < 1.{
            start_points_mu_ik.extend( linspace(0.01, 2., 50).step_by(1));
        } else if mean_mu < 10. && min == 0. {
            start_points_mu_ik.extend( linspace(min, max, 50).step_by(1));
        } else if min == 0. {
            min = mean_mu - mean_mu * 0.1;
            start_points_mu_ik.extend( linspace(min, max, 50).step_by(1));
        } else {
            start_points_mu_ik.extend( linspace(min, max, 50).step_by(1));
        }
        start_points_mu_ik.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // round each value to 3 decimal places
        start_points_mu_ik = start_points_mu_ik.iter().map(|x| (x * 1000.).round() / 1000.).collect();
        start_points_mu_ik.dedup();
        // start_points_mu_ik = start_points_mu_ik.iter().step_by(2).map(|x| *x).collect();


        let mut possible_f: Vec<f64> = linspace(0.05, 5., 100).step_by(1).collect();
        possible_f.extend( linspace(5., 10., 20).step_by(1));
        possible_f.extend( linspace(10.5, 20., 20).step_by(1));
        // println!("len possible_f {:?}", possible_f.len());
        possible_f.sort_by(|a, b| a.partial_cmp(b).unwrap());
        possible_f.dedup();


        let mut thetas: Vec<f64> = linspace(0.01, 0.1, 10).take(1).collect();
        thetas.extend( linspace(0.1, 1., 10).step_by(1));
        thetas.extend( linspace(1.5, 10., 18).step_by(2));
        thetas.extend( linspace(11., 165., 155).step_by(10));
        // println!("len thetas {:?}", thetas.len());
        thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        thetas.dedup();

        let mut all_mu_ik =  start_points_mu_ik.clone();
        for (f, mu_1) in iproduct!(possible_f.clone(), start_points_mu_ik.clone()) {
            let mut mu_2 = f * (mu_1 + c) - c;
            if mu_2 > 0. {
                if mu_2 < 0.1 {
                    // round mu_2 to 3 decimals
                    mu_2 = (mu_2 * 1000.).round() / 1000.;
                } else if mu_2 < 100. {
                    // round mu_2 to 2 decimals
                    mu_2 = (mu_2 * 100.).round() / 100.;
                } else {
                     // round mu_2 to 1 decimal
                    mu_2 = (mu_2 * 10.).round() / 10.;
                }
                all_mu_ik.push(mu_2);
            }
        }
        all_mu_ik.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
        all_mu_ik.dedup();
        // println!("len all_mu_ik {:?}", all_mu_ik.len());

        Ok(QueryPoints {
            start_points_mu_ik,
            all_mu_ik,
            thetas,
            possible_f,
        })
    }

}

fn compute_min_max_values(mean_disp_estimates: &HashMap<String, Estimates>) -> HashMap<usize, (f64, f64)> {
    // Initialize the HashMap to store the minimum and maximum values for each feature
    let mut min_max_values: HashMap<usize, (f64, f64)> = HashMap::new();

    // Iterate over each sample
    for (_, estimates) in mean_disp_estimates {
        // Iterate over each feature in the estimate
        for (i, &value) in estimates.means().iter().enumerate() {
            // Update the minimum and maximum values for the current feature
            let (min_value, max_value) = min_max_values.entry(i).or_insert((value, value));
            if value < *min_value {
                *min_value = value;
            }
            if value > *max_value {
                *max_value = value;
            }
        }
    }
    min_max_values
}


// calc query points for each feature
pub(crate) fn calc_query_points(c: f64, mean_disp_estimates: HashMap<String, Estimates> , sample_ids: Vec<String>, feature_ids: Array1<String>, id: usize ) ->QueryPoints {
    // pub(crate) fn calc_query_points(c: f64, mean_disp_estimates: HashMap<String, Estimates> , sample_ids: Vec<String>, feature_ids: Array1<String> ) ->HashMap<String, QueryPoints> {
//     // let mut query_points_per_feature = HashMap::<String, QueryPoints>::new();
   // get minimum, mean and maximum mu_ik per feature

   let mut means_per_feature : Array1<f64> = Array1::zeros(Dim([mean_disp_estimates[&sample_ids[0]].means().len()]));
   means_per_feature = sample_ids
       .iter()
       .fold(means_per_feature, |acc: Array1<f64>, x: &String| acc + mean_disp_estimates[x].means());

    let min_max_values = compute_min_max_values(&mean_disp_estimates);
    let number_of_samples = sample_ids.len() as f64;
    means_per_feature = means_per_feature.map(|x| -> f64 {x / number_of_samples});


//    let mut min_per_feature : Array1<f64> = Array1::zeros(Dim([mean_disp_estimates[&sample_ids[0]].means().len()]));
//    // get list of minimum mu_ik per feature
//     min_per_feature =  sample_ids
//        .iter()
//        .fold(min_per_feature, |acc: Array1<f64>, x: &String| {
//            let mut min = acc.clone();
//            for (i, x) in mean_disp_estimates[x].means().iter().enumerate() {
//                if x.is_finite() {
//                    if x < &min[i] {
//                        min[i] = *x;
//                    }
//                }
//            }
//            min
//        });
    // println!("min_per_feature 1 {:?}", min_per_feature[190432]);
    // println!("min_per_feature -1 {:?}", min_per_feature[min_per_feature.len()-1]);
    //enumerate over features
    // for (i, feature_id ) in feature_ids.iter().enumerate() {
    //     let mean = means_per_feature[i];
    //     let query_points = QueryPoints::new(c, mean, *min_max_values.get(&i).unwrap()).unwrap();
    //     query_points_per_feature.insert(feature_id.to_string(), query_points);
    // }

    // query_points_per_feature
    let mean = means_per_feature[id];
    let query_points = QueryPoints::new(c, mean, *min_max_values.get(&id).unwrap()).unwrap();
    query_points

}
