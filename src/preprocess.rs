use std::collections::HashMap;
use std::fs::File;
use std::io::stdout;
use std::path::{Path, PathBuf};

use anyhow::Result;
use getset::Getters;
use ndarray::{s, Array1, Axis};
use ndarray_stats::{interpolate, Quantile1dExt, QuantileExt};
use noisy_float::types::N64;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};

use crate::errors::Error;
use crate::kallisto::{self, KallistoQuant};

pub(crate) fn preprocess(kallisto_quants: &[PathBuf], sample_ids: &[String]) -> Result<()> {
    if kallisto_quants.len() < 2 {
        return Err(Error::NotEnoughQuants.into());
    }

    let quants: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|kallisto_quant| KallistoQuant::new(kallisto_quant))
        .collect();
    let quants = quants?;

    let scale_factors = scale_factors(&quants, sample_ids)?;
    dbg!(&scale_factors);

    let mean_disp_estimates = mean_disp_estimates(&quants, sample_ids)?;
    let preprocessing = Preprocessing {
        scale_factors,
        mean_disp_estimates,
        feature_ids: quants[0].feature_ids()?,
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
}

impl Preprocessing {
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(Preprocessing::deserialize(&mut Deserializer::new(
            File::open(path)?,
        ))?)
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

fn scale_factors(
    kallisto_quants: &[KallistoQuant],
    sample_ids: &[String],
) -> Result<HashMap<String, f64>> {
    // TODO do we really need to norm counts by the len?

    let counts: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|quant| quant.len_norm_counts())
        .collect();
    let mut counts = counts?;

    let upper_quartiles: Array1<N64> = counts
        .iter_mut()
        .map(|feature_counts| {
            feature_counts
                .quantile_mut(N64::unchecked_new(0.75), &interpolate::Linear)
                .unwrap()
                .clone()
        })
        .collect();

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

#[derive(Serialize, Deserialize, Debug, Getters)]
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
