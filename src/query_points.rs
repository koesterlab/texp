use std::collections::HashMap;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use anyhow::Result;
use getset::Getters;
use itertools::iproduct;
use itertools_num::linspace;
use serde_derive::{Deserialize, Serialize};

use crate::preprocess::Estimates;
use crate::preprocess::Preprocessing;
use crate::prior::{Prior, PriorParameters};

#[derive(Serialize, Deserialize, Debug, Getters)]
pub(crate) struct QueryPoints {
    #[get = "pub(crate)"]
    start_points_mu_ik: Vec<f64>,
    #[get = "pub(crate)"]
    all_mu_ik: Vec<f64>,
    #[get = "pub(crate)"]
    thetas: Vec<f64>,
    #[get = "pub(crate)"]
    possible_f: Vec<f64>,
}

#[derive(Debug)]
pub(crate) struct QueryPointsPerFeature {
    per_feature: Vec<QueryPoints>,
}

impl QueryPointsPerFeature {
    #[inline]
    pub(crate) fn get(&self, feature_id: usize) -> &QueryPoints {
        &self.per_feature[feature_id]
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.per_feature.len()
    }
}

impl QueryPointsPerFeature {
    pub(crate) fn new(preprocessing: &Preprocessing, c: f64) -> Self {
        let sample_ids: Vec<_> = preprocessing.scale_factors().keys().cloned().collect();

        let num_features = preprocessing.feature_ids().len();
        let num_samples = sample_ids.len() as f64;

        // ---- mean per feature ----
        let mut means_per_feature = vec![0.0; num_features];

        for sample_id in &sample_ids {
            let est = &preprocessing.mean_disp_estimates()[sample_id];
            for (i, &v) in est.means().iter().enumerate() {
                means_per_feature[i] += v;
            }
        }

        for v in &mut means_per_feature {
            *v /= num_samples;
        }


        // ---- standard deviation per feature ----
        let mut std_per_feature = vec![0.0; num_features];

        for sample_id in &sample_ids {
            let est = &preprocessing.mean_disp_estimates()[sample_id];
            for (i, &v) in est.means().iter().enumerate() {
                let diff = v - means_per_feature[i];
                std_per_feature[i] += diff * diff;
            }
        }

        for std in &mut std_per_feature {
            *std = (*std / num_samples).sqrt();
        }

        // ---- min / max per feature ----
        let min_max_values = compute_min_max_values(preprocessing.mean_disp_estimates());

        // // ---- build QueryPoints for all features ----
        // let per_feature = (0..num_features)
        //     .map(|i| QueryPoints::new(c, means_per_feature[i], min_max_values[&i]).unwrap())
        //     .collect();

         // ---- build QueryPoints using mean ± std ----
        let per_feature = (0..num_features)
            .map(|i| {
                let min = (means_per_feature[i] - std_per_feature[i]).max(0.0);
                let max = (means_per_feature[i] + std_per_feature[i]).min(min_max_values[&i].1);

                QueryPoints::new(preprocessing, c, means_per_feature[i], (min, max)).unwrap()
            })
            .collect();

        Self { per_feature }
    }
}

impl QueryPoints {
    // get query points for one feature
    pub(crate) fn new(preprocessing: &Preprocessing, c: f64, mean_mu: f64, min_max: (f64, f64)) -> Result<Self> {
        let mut min = min_max.0;
        let max = min_max.1;
        let mut start_points_mu_ik = vec![0.];
        // if mean_mu < 0.1 {
        //     start_points_mu_ik.extend(linspace(0.001, 0.2, 40).step_by(1));
        // } else if mean_mu < 1. {
        //     start_points_mu_ik.extend(linspace(0.01, 2., 40).step_by(1));
        // } else if mean_mu < 10. && min == 0. {
        //     start_points_mu_ik.extend(linspace(min, max, 40).step_by(1));
        // } else if min == 0. {
        //     min = mean_mu - mean_mu * 0.1;
        //     start_points_mu_ik.extend(linspace(min, max, 40).step_by(1));
        // } else {
        //     start_points_mu_ik.extend(linspace(min, max, 40).step_by(1));
        // }
        start_points_mu_ik.extend(linspace(min, max, 20).step_by(1));
        start_points_mu_ik.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // round each value to 3 decimal places
        start_points_mu_ik = start_points_mu_ik
            .iter()
            .map(|x| (x * 100.).round() / 100.)
            .collect();
        start_points_mu_ik.dedup();
        // start_points_mu_ik = start_points_mu_ik.iter().step_by(2).map(|x| *x).collect();

        let mut possible_f: Vec<f64> = linspace(0.05, 4., 25).step_by(1).collect();
        possible_f.extend(linspace(4., 10., 11).step_by(1));
        // possible_f.extend(linspace(10.5, 20., 10).step_by(1));
        // println!("len possible_f {:?}", possible_f.len());
        possible_f.sort_by(|a, b| a.partial_cmp(b).unwrap());
        possible_f.dedup();

        let prior = Prior::new(&preprocessing.prior_parameters())?;
        // fixed seed
        let mut rng = StdRng::seed_from_u64(100);

        //  draw 20 dispersions
        let mut thetas_rand = prior.sample_n(20, &mut rng);
        thetas_rand.extend([0.05, 0.1, 0.25, 0.5]);
        thetas_rand.sort_by(|a, b| a.partial_cmp(b).unwrap());
        thetas_rand.dedup();
        // dbg!(&thetas_rand);
        // dbg!(thetas_rand.len());

        // let mut thetas: Vec<f64> = linspace(0.1, 0.1, 1).collect();
        // thetas.extend(linspace(0.1, 1., 10).step_by(1));
        // thetas.extend(linspace(1.5, 10., 15).step_by(2));
        // thetas.extend(linspace(11., 150., 115).step_by(10));
        // // println!("len thetas {:?}", thetas.len());
        // thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // thetas.dedup();
        // println!("len thetas {:?}", thetas.len());
        // println!("len thetas_rand {:?}", thetas_rand.len());
        let thetas = thetas_rand;

        let mut all_mu_ik = start_points_mu_ik.clone();
        for (f, mu_1) in iproduct!(possible_f.clone(), start_points_mu_ik.clone()) {
            let mut mu_2 = f * (mu_1 + c) - c;
            if mu_2 > 0. {
                if mu_2 < 0.1 {
                    // round mu_2 to 3 decimals
                    mu_2 = (mu_2 * 100.).round() / 100.;
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

fn compute_min_max_values(
    mean_disp_estimates: &HashMap<String, Estimates>,
) -> HashMap<usize, (f64, f64)> {
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

fn compute_means_per_feature(preprocessing: &Preprocessing) -> Vec<f64> {
    let sample_ids: Vec<_> = preprocessing.scale_factors().keys().cloned().collect();
    let num_features = preprocessing.feature_ids().len();
    let num_samples = sample_ids.len() as f64;

    let mut means = vec![0.0; num_features];

    for sample_id in sample_ids {
        let est = &preprocessing.mean_disp_estimates()[&sample_id];
        for (i, &v) in est.means().iter().enumerate() {
            means[i] += v;
        }
    }

    for v in &mut means {
        *v /= num_samples;
    }

    means
}
