
use std::collections::HashMap;
use std::fs::File;
use std::io::stdout;
use std::path::{Path, PathBuf};

use anyhow::Result;
use getset::Getters;
use ndarray::{Array1, Axis, Dim};
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};

use crate::preprocess::Estimates;


#[derive(Serialize, Deserialize, Debug, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct GroupMeans {
    group_means: Array1<f64>,
}

impl GroupMeans {
    pub(crate) fn from_path(path: &Path) -> Result<Self> {
        Ok(GroupMeans::deserialize(&mut Deserializer::new(
            File::open(path)?,
        ))?)
    }
}

// let group_means = calc_group_means(&mean_disp_estimates, sample_ids)?;


pub(crate) fn calc_group_means(
    mean_disp_estimates: &HashMap<String, Estimates>,
    sample_ids: &[String],
) -> Result<()> {
    dbg!("function group_means");

    let mut means : Array1<f64> = Array1::zeros(Dim([mean_disp_estimates[&sample_ids[0]].means().len()]));
    means = sample_ids
        .iter()
        .fold(means, |acc: Array1<f64>, x: &String| acc + mean_disp_estimates[x].means());
    let number_of_samples = sample_ids.len() as f64;
    means = means.map(|x| -> f64 {x / number_of_samples});
    // Ok(means)


    let group_means = GroupMeans{group_means : means};
    group_means.serialize(&mut Serializer::new(stdout()))?;

    Ok(())


}