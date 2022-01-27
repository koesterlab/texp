use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N32;
use rayon::prelude::*;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;

use crate::common::{window, Outdir, ProbDistribution};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prior::Prior;
use crate::sample_expression::SampleInfo;
// use crate::group_expression::group_expression


pub(crate) fn diff_exp(
    c: usize,
//     preprocessing: &Path,
    group_expression_paths: &[PathBuf],
    out_dir: &Path,
) -> Result<()> {

    println!("{}", c);
    Ok(())

}
