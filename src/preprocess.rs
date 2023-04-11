//! This infers scale factors, mean and dispersion from Kallisto results.

use std::collections::HashMap;
use std::fs::File;
use std::io::stdout;
use std::path::{Path, PathBuf};
use std::thread;

use anyhow::Result;
use getset::Getters;
use ndarray::{Array1, Axis, Dim};
use ndarray_stats::{interpolate, Quantile1dExt, QuantileExt};
use noisy_float::types::N64;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};
use itertools_num::linspace;
use itertools::iproduct;

use crate::errors::Error;
use crate::kallisto::KallistoQuant;
use crate::prior::{Prior, PriorParameters};

pub(crate) fn preprocess(
    c :f64,
    kallisto_quants: &[PathBuf],
    sample_ids: &[String],
    prior_parameters: PriorParameters,    
) -> Result<()> {
    if kallisto_quants.len() < 1 {
        return Err(Error::NotEnoughQuants.into());
    }
    let quants: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|kallisto_quant| KallistoQuant::new(kallisto_quant))
        .collect();
    let quants = quants?;

    let mut scale_factors: HashMap<String, f64>;
    if kallisto_quants.len() > 1 {
        scale_factors = calc_scale_factors(&quants, sample_ids)?;
    } else {
        scale_factors =  HashMap::new();
        scale_factors.insert((*sample_ids[0]).to_string(), 1.);
    }
    dbg!(&scale_factors);

    let mean_disp_estimates = mean_disp_estimates(&quants, sample_ids)?;

    let feature_ids = quants[0].feature_ids()?;

    let query_points = calc_query_points(c, mean_disp_estimates.clone(), sample_ids, feature_ids.clone());

    let preprocessing = Preprocessing {
        scale_factors,
        mean_disp_estimates,
        query_points,
        feature_ids: feature_ids,
        prior_parameters,
    };

    preprocessing.serialize(&mut Serializer::new(stdout()))?;

    Ok(())
}

#[derive(Serialize, Deserialize, Debug, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct Preprocessing {
    scale_factors: HashMap<String, f64>,
    mean_disp_estimates: HashMap<String, Estimates>,
    query_points: HashMap<String, QueryPoints>,
    feature_ids: Array1<String>,
    prior_parameters: PriorParameters,
}

impl Preprocessing {
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(Preprocessing::deserialize(&mut Deserializer::new(
            File::open(path)?,
        ))?)
    }

    pub(crate) fn prior(&self) -> Result<Prior> {
        Prior::new(self.prior_parameters())
    }

    // pub(crate) fn prior(&self) -> Result<Prior> {
    //     Prior::new(0.25)
    // }

    pub(crate) fn interpolate_dispersion(&self, feature_idx: usize) -> Option<f64> {
        let disp = |estimates: &Estimates| estimates.dispersions[feature_idx];
        let count = self.mean_disp_estimates.values().filter_map(&disp).count();
        if count == 0 {
            None
        } else {
            Some(
                self.mean_disp_estimates
                    .values()
                    .filter_map(&disp)
                    .sum::<f64>()
                    / count as f64,
            )
        }
    }
}

fn calc_scale_factors(
    kallisto_quants: &[KallistoQuant],
    sample_ids: &[String],
) -> Result<HashMap<String, f64>> {
    // TODO do we really need to norm counts by the len?

    let counts: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|quant| quant.len_norm_counts())
        .collect();
    let mut counts = counts?;

    let child = thread::Builder::new().stack_size(32 * 1024 * 1024).spawn(move || {
        let upper_quartiles: Array1<N64> = counts
        .iter_mut()
        .map(|feature_counts| {
            feature_counts
                .quantile_mut(N64::unchecked_new(0.75), &interpolate::Linear)
                .unwrap()
                .clone()
        })
        .collect();
        return upper_quartiles;
    }).unwrap();

    let upper_quartiles = child.join().unwrap();

    // let upper_quartiles: Array1<N64> = counts
    //     .iter_mut()
    //     .map(|feature_counts| {
    //         feature_counts
    //             .quantile_mut(N64::unchecked_new(0.75), &interpolate::Linear)
    //             .unwrap()
    //             .clone()
    //     })
    //     .collect();
    let max_quartile = upper_quartiles.max()?.clone();
    let scale_factors = upper_quartiles.mapv(|quartile| max_quartile / quartile);
    Ok(sample_ids
        .iter()
        .cloned()
        .zip(
            scale_factors
                .iter()
                .map(|scale_factor| (*scale_factor).into()),
        )
        .collect())
    
}

#[derive(Serialize, Deserialize, Debug, Getters, Clone)]
#[getset(get = "pub(crate)")]
pub(crate) struct Estimates {
    dispersions: Array1<Option<f64>>,
    means: Array1<f64>,
}


impl Estimates {
    fn new(kallisto_quant: &KallistoQuant) -> Result<Self> {
        let bootstrapped_counts = kallisto_quant.bootstrapped_counts()?;
        let means = bootstrapped_counts.mean_axis(Axis(0)).unwrap();
        let stds = bootstrapped_counts.std_axis(Axis(0), 1.0);
        let dispersions = stds / &means;

        Ok(Estimates {
            dispersions: dispersions.mapv(|d| if d.is_nan() { None } else { Some(d) }),
            means,
        })
    }
}

fn mean_disp_estimates(
    kallisto_quants: &[KallistoQuant],
    sample_ids: &[String],
) -> Result<HashMap<String, Estimates>> {
    let estimates: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|quant| Estimates::new(quant))
        .collect();

    Ok(sample_ids.iter().cloned().zip(estimates?).collect())
}


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
        // let mut start_points_mu_ik: Vec<f64> = linspace(0., 500., 50).take(1).collect();
        
        // } else if mean_mu < 10. {
        //     start_points_mu_ik.extend( linspace(1., 10., 50).step_by(1));
        // } else if mean_mu < 100. {
        //     start_points_mu_ik.extend( linspace(10., 100., 50).step_by(1));
        // } else  if mean_mu < 500. {
        //     start_points_mu_ik.extend( linspace(100., 500., 50).step_by(1));
        // } else  if mean_mu < 1000. {
        //     start_points_mu_ik.extend( linspace(500., 1000., 50).step_by(1));
        // } else  if mean_mu < 5000. {
        //     start_points_mu_ik.extend( linspace(1000., 5000., 50).step_by(1));
        // } else  if mean_mu < 10000. {
        //     start_points_mu_ik.extend( linspace(5000., 10000., 50).step_by(1));
        // } else {
        //     start_points_mu_ik.extend( linspace(10000., 200 000., 50).step_by(1));
        // }
        // let mut start_points_mu_ik: Vec<f64> = linspace(0., 1., 151).collect();
        // start_points_mu_ik.extend( linspace(1.5, 100., 198).step_by(10));
        // start_points_mu_ik.extend( linspace(101., 500., 200).step_by(20));
        // start_points_mu_ik.extend( linspace(600., 3000., 25));
        // start_points_mu_ik.extend( linspace(4000., 10000., 7));
        // println!("len start_points_mu_ik {:?}", start_points_mu_ik.len());
        start_points_mu_ik.sort_by(|a, b| a.partial_cmp(b).unwrap());
        start_points_mu_ik.dedup();


        let mut possible_f: Vec<f64> = linspace(0.05, 5., 100).step_by(1).collect();
        possible_f.extend( linspace(5., 10., 40).step_by(1));
        possible_f.extend( linspace(10.5, 20., 30).step_by(1));
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
            let mu_2 = f * (mu_1 + c) - c;
            if mu_2 > 0. {
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
fn calc_query_points(c: f64, mean_disp_estimates: HashMap<String, Estimates> , sample_ids: &[String], feature_ids: Array1<String> ) ->HashMap<String, QueryPoints> {
    let mut query_points_per_feature = HashMap::<String, QueryPoints>::new();
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
    for (i, feature_id ) in feature_ids.iter().enumerate().skip(190432) {
        let mean = means_per_feature[i];
        let query_points = QueryPoints::new(c, mean, *min_max_values.get(&i).unwrap()).unwrap();
        query_points_per_feature.insert(feature_id.to_string(), query_points);
    }

    query_points_per_feature

}
