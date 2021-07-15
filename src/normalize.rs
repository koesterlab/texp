use std::io::stdout;
use std::path::PathBuf;

use anyhow::Result;
use ndarray::Array1;
use ndarray_stats::{interpolate, Quantile1dExt, QuantileExt};
use noisy_float::types::N64;
use rmp_serde::Serializer;
use serde::Serialize;
use serde_derive::Serialize;

use crate::errors::Error;
use crate::kallisto::KallistoQuant;

pub(crate) fn normalize(kallisto_quants: &[PathBuf], sample_ids: &[String]) -> Result<()> {
    if kallisto_quants.len() < 2 {
        return Err(Error::NotEnoughQuants.into());
    }

    // TODO do we really need to norm counts by the len?
    let counts: Result<Vec<_>> = kallisto_quants
        .iter()
        .map(|kallisto_quant| -> Result<Array1<N64>> {
            let quant = KallistoQuant::new(kallisto_quant)?;
            quant.len_norm_counts()
        })
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

    let scale_factors: Vec<_> = scale_factors
        .iter()
        .zip(sample_ids)
        .map(|(scale_factor, sample_id)| ScaleFactor {
            sample_id: sample_id.clone(),
            scale_factor: scale_factor.clone().into(),
        })
        .collect();

    scale_factors.serialize(&mut Serializer::new(stdout()))?;

    Ok(())
}

#[derive(Serialize)]
struct ScaleFactor {
    sample_id: String,
    scale_factor: f64,
}
