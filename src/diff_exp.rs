use std::fs;
use std::path::{Path, PathBuf};
use std::mem;

use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N64;
use rayon::prelude::*;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use itertools_num::linspace;


use crate::common::{window, Outdir, ProbDistribution};
use crate::errors::Error;
use crate::preprocess::Preprocessing;
use crate::prior::Prior;
use crate::sample_expression::SampleInfo;
// use crate::group_expression::group_expression


pub(crate) fn diff_exp(
    c: f64,
//     preprocessing: &Path,
    group_expression_paths: &[PathBuf],
    out_dir: &Path,
) -> Result<()> {

    let mut prob_dist_i_k1 = ProbDistribution::default();
    prob_dist_i_k1.insert(N64::new(1. as f64), LogProb(10.0f64.ln()));
    let mut prob_dist_i_k2 = ProbDistribution::default();
    prob_dist_i_k2.insert(N64::new(1. as f64), LogProb(10.0f64.ln()));

    let mut max_pos_i_k1 = prob_dist_i_k1.get_max_position().unwrap();
    let mut max_pos_i_k2 = prob_dist_i_k2.get_max_position().unwrap();

    let range_for_f = vec![1./5., 1./4., 1./3., 1./2., 1., 2., 3., 4., 5.];
    for i in &range_for_f {
        println!("{}", i);
    }
    if max_pos_i_k2 < max_pos_i_k1 {
        mem::swap(&mut prob_dist_i_k1, &mut prob_dist_i_k2);
        mem::swap(&mut max_pos_i_k1, &mut max_pos_i_k2);
    }


    let mut diff_exp_distribution = ProbDistribution::default();

    range_for_f
        .par_iter()
        .try_for_each(|f|-> Result<()> {
            let f_max_i_k2 = (-f * c + c + max_pos_i_k1.raw()) / f;
            let range_for_x = linspace(0., max_pos_i_k1.raw(), 15).chain(linspace(max_pos_i_k1.raw(), f_max_i_k2, 15)).chain(linspace(f_max_i_k2, f_max_i_k2 + 15., 15));

            let density = |i, x| {
                prob_dist_i_k1.get(x) + prob_dist_i_k2.get(&N64::new(f * (x.raw() + c) -c));
            };

            let mut value = LogProb::ln_simpsons_integrate_exp(
                density,
                0.f64,
                15.f64,
                11,
            );
            diff_exp_distribution.insert(*f, value);
           Ok(())
       });

   

    Ok(())

