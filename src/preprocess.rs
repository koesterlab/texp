//! This infers scale factors, mean and dispersion from Kallisto results.

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::collections::HashMap;
use std::fs::File;
use std::io::stdout;
use std::path::{Path, PathBuf};
use std::thread;

use anyhow::Result;
use getset::Getters;
use itertools_num::linspace;
use ndarray::{Array1, Axis};
use ndarray_stats::{interpolate, Quantile1dExt, QuantileExt};
use noisy_float::types::N64;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};
use statrs::function::beta::ln_beta;

use crate::errors::Error;
use crate::kallisto::KallistoQuant;
use crate::prior::{Prior, PriorParameters};

pub(crate) fn preprocess(
    _c: f64,
    kallisto_quants: &[PathBuf],
    sample_ids: &[String],
    prior_parameters: PriorParameters,
) -> Result<()> {
    if kallisto_quants.is_empty() {
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
        scale_factors = HashMap::new();
        scale_factors.insert((*sample_ids[0]).to_string(), 1.);
    }
    dbg!(&scale_factors);

    let mean_disp_estimates = mean_disp_estimates(&quants, sample_ids)?;

    let feature_ids = quants[0].feature_ids()?;

    //TODO Get thetas from query points instead of hardcoding them here. We need to make sure that the same thetas are used for both the cache and the query points.
    let mut thetas: Vec<f64> = linspace(0.01, 0.1, 5).collect();
    thetas.extend(linspace(0.1, 1., 10).step_by(1));
    thetas.extend(linspace(1.5, 10., 15).step_by(2));
    thetas.extend(linspace(11., 165., 115).step_by(10));
    // println!("len thetas {:?}", thetas.len());
    thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
    thetas.dedup();
    dbg!(&thetas);
    dbg!(thetas.len());

    let ln_beta_caches = thetas
        .iter()
        .map(|&theta| LnBetaCache::new(theta, 10000))
        .collect();

    let prior = Prior::new(&prior_parameters)?;
    // fixed seed
    let mut rng = StdRng::seed_from_u64(12345);

    // draw 150 dispersions
    let mut thetas_rand = prior.sample_n(35, &mut rng);
    thetas_rand.sort_by(|a, b| a.partial_cmp(b).unwrap());
    thetas_rand.dedup();
    dbg!(&thetas_rand);
    dbg!(thetas_rand.len());

    let preprocessing = Preprocessing {
        scale_factors,
        mean_disp_estimates,
        feature_ids,
        prior_parameters,
        ln_beta_caches,
    };

    preprocessing.serialize(&mut Serializer::new(stdout()))?;

    Ok(())
}

#[derive(Serialize, Deserialize, Debug, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct Preprocessing {
    scale_factors: HashMap<String, f64>,
    mean_disp_estimates: HashMap<String, Estimates>,
    feature_ids: Array1<String>,
    prior_parameters: PriorParameters,
    ln_beta_caches: Vec<LnBetaCache>, // one per theta
}

impl Preprocessing {
    //constructor for preprocessing
    pub(crate) fn new(
        scale_factors: HashMap<String, f64>,
        mean_disp_estimates: HashMap<String, Estimates>,
        feature_ids: Array1<String>,
        prior_parameters: PriorParameters,
        ln_beta_caches: Vec<LnBetaCache>,
    ) -> Self {
        Preprocessing {
            scale_factors,
            mean_disp_estimates,
            feature_ids,
            prior_parameters,
            ln_beta_caches,
        }
    }

    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(Preprocessing::deserialize(&mut Deserializer::new(
            File::open(path)?,
        ))?)
    }

    pub(crate) fn prior(&self) -> Result<Prior> {
        Prior::new(self.prior_parameters())
    }

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

    let child = thread::Builder::new()
        .stack_size(32 * 1024 * 1024)
        .spawn(move || {
            let upper_quartiles: Array1<N64> = counts
                .iter_mut()
                .map(|feature_counts| {
                    feature_counts
                        .quantile_mut(N64::unchecked_new(0.75), &interpolate::Linear)
                        .unwrap()
                })
                .collect();
            upper_quartiles
        })
        .unwrap();

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
    let max_quartile = *upper_quartiles.max()?;
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

    pub(crate) fn new_from_arrays(dispersions: Array1<Option<f64>>, means: Array1<f64>) -> Self {
        Estimates { dispersions, means }
    }
}

fn mean_disp_estimates(
    kallisto_quants: &[KallistoQuant],
    sample_ids: &[String],
) -> Result<HashMap<String, Estimates>> {
    let estimates: Result<Vec<_>> = kallisto_quants.iter().map(Estimates::new).collect();

    Ok(sample_ids.iter().cloned().zip(estimates?).collect())
}

#[derive(Serialize, Deserialize, Debug, Getters, Clone)]
pub struct LnBetaCache {
    n: f64,
    values: Vec<f64>, // ln_beta(x+1, n)
}

impl LnBetaCache {
    pub fn new(theta: f64, initial_x: usize) -> Self {
        let n = 1.0 / theta;
        let mut v = Vec::with_capacity(initial_x + 1);
        for x in 0..=initial_x {
            v.push(ln_beta((x as f64) + 1.0, n));
        }
        Self { n, values: v }
    }

    #[inline]
    pub fn get(&self, x: usize) -> f64 {
        self.values[x]
    }
}
