//! reduce the number of features in a dataset to only those specified in a list

use anyhow::Result;
use ndarray::{Array, Array1, s};
use std::path::Path;
use std::collections::HashMap;
use serde::Serialize as SerdeSerialize;
use rmp_serde::{Deserializer, Serializer};
use std::io::stdout;

// use crate::common::Outdir; 
use crate::preprocess::{Preprocessing, Estimates};

pub(crate) fn reduce_features(
    preprocessing: &Path,
    feature_ids: &Path,
    // out_dir_path: &Path,
) -> Result<()> {

    // let out_dir = Outdir::create(out_dir_path)?;

    let preprocessing = Preprocessing::from_path(preprocessing)?;
    // read list of feature ids in vector
    let wanted_feature_ids = std::fs::read_to_string(feature_ids)?;
    // skip empty strings
    let wanted_feature_ids: Vec<&str> = wanted_feature_ids.split("\n").collect();
    let wanted_feature_ids: Vec<&str> = wanted_feature_ids.iter().filter(|&s| !s.is_empty()).map(|&s| s).collect();
    // println!("wanted_feature_ids {:?}", wanted_feature_ids);
    // let reduced_preprocessing = preprocessing.reduce_features(wanted_feature_ids)?;
    let scale_factors = preprocessing.scale_factors();
    // mean_disp_estimates = Hasmap of sample_id -> esitmates
    let mean_disp_estimates = preprocessing.mean_disp_estimates();
    let prior_parameters = preprocessing.prior_parameters();
    let feature_ids = preprocessing.feature_ids();

    //get indices of wanted features in feature_ids
    let mut wanted_feature_indices: Vec<usize> = Vec::new();
    for wanted_feature_id in &wanted_feature_ids {
        let wanted_feature_index = feature_ids.iter().position(|x| *x == *wanted_feature_id);
        match wanted_feature_index {
            Some(i) => wanted_feature_indices.push(i),
            None => println!("feature {} not found", wanted_feature_id),
        }
    }

    // filter out unwanted features from mean_disp_estimates for each sample
    let mut mean_disp_estimates_filtered: HashMap<String, Estimates> = HashMap::new();
    for (sample_id, estimates) in mean_disp_estimates.iter() {
        let mut dispersions = Vec::new();
        let mut means = Vec::new();
        
        for wanted_feature_index in wanted_feature_indices.iter() {
            dispersions.push(estimates.dispersions()[*wanted_feature_index]);
            means.push(estimates.means()[*wanted_feature_index]);
        }
        let dispersions = Array1::from(dispersions);
        let means = Array1::from(means);
        let estimates_filtered = Estimates::new_from_arrays(dispersions, means);
        mean_disp_estimates_filtered.insert(sample_id.to_string(), estimates_filtered);
    }


    //get array of string from wanted_feature_ids str vec
    let feature_ids: Array1<String> = wanted_feature_ids
        .iter()
        .map(|&s| s.to_string())
        .collect::<Vec<String>>()
        .into();

    let filtered_preprocessing = Preprocessing::new(
        scale_factors.clone(),
        mean_disp_estimates_filtered,
        feature_ids,
        *prior_parameters
    );
    filtered_preprocessing.serialize(&mut Serializer::new(stdout()))?;

    Ok(())
}