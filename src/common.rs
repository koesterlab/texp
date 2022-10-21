use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use derefable::Derefable;
use itertools_num::linspace;
// use noisy_float::prelude::Float;
// use noisy_float::types::N32;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;

use crate::errors::Error;


pub(crate) struct Pair {
    pub left: f64,
    pub right: f64,
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct Point {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug)]
pub(crate) struct Square {
    pub top_left: usize,
    pub top_right: usize,
    pub bot_left: usize,
    pub bot_right: usize,
}

/// Return left and right window around given parameter mean
#[allow(unused)]
pub(crate) fn window(mean: f64) -> (impl Iterator<Item = f64>, impl Iterator<Item = f64>) {
    // TODO: think about larger steps, binary search etc. to optimize instead of just having a fixed number of steps.
    (
        linspace(mean / 5.0, mean, 20).rev().skip(1),
        linspace(mean, 10000.0 * mean, 120),
    )
}

#[allow(unused)]
pub(crate) fn window_x(mean: f64) -> (impl Iterator<Item = f64>, impl Iterator<Item = f64>) {
    // TODO: think about larger steps, binary search etc. to optimize instead of just having a fixed number of steps.
    (linspace(0.0, 50., 51), linspace(51., 500., 450))
}

// pub(crate) fn window_f(
//     max_prob_fold_change: f64,
// ) -> (impl Iterator<Item = f64>, impl Iterator<Item = f64>) {
//     // TODO: think about larger steps, binary search etc. to optimize instead of just having a fixed number of steps.
//     (
//         linspace(max_prob_fold_change / 5.0, max_prob_fold_change, 20)
//             .rev()
//             .skip(1),
//         linspace(max_prob_fold_change, 30.0 * max_prob_fold_change, 30),
//     )
// }

pub(crate) fn difference_to_big( // If difference is > than 1/10.000 of maximum probability
    estimated_value: LogProb,
    calculated_value: LogProb,
    prob_dist_max: LogProb
) -> bool {
    // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
    if estimated_value > calculated_value {
        // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
        // println!("ln_sub-max1 {:?}, thresh {:?}, diffToBig {:?}", estimated_value.ln_sub_exp(calculated_value)- prob_dist_max, LogProb::from(0.01_f64.ln()),estimated_value.ln_sub_exp(calculated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()) );
        return estimated_value.ln_sub_exp(calculated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()); 
    }else {
        // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
        // println!("ln_sub-max2 {:?}, thresh {:?}, diffToBig {:?}", calculated_value.ln_sub_exp(estimated_value) - prob_dist_max, LogProb::from(0.01_f64.ln()), calculated_value.ln_sub_exp(estimated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()) );
        return calculated_value.ln_sub_exp(estimated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()); 
    }
    // if (estimated_value.exp() - calculated_value.exp()).abs() / prob_dist_max > 0.0001 
}

// pub(crate) fn interpolate_pmf(
//     value: N32,
//     lower: N32,
//     upper: N32,
//     prob_lower: LogProb,
//     prob_upper: LogProb,
// ) -> LogProb {
//     let len = upper - lower;
//     (prob_lower + LogProb(f64::from(((upper - value) / len).ln())))
//         .ln_add_exp(prob_upper + LogProb(f64::from(((value - lower) / len).ln())))
// }

//--------------------OutDir--------------------

#[derive(Derefable)]
pub(crate) struct Outdir {
    #[deref]
    path: PathBuf,
}

impl Outdir {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        if !path.exists() {
            return Err(Error::NotExistingOutputDir {
                path: path.to_owned(),
            }
            .into());
        }
        Ok(Outdir {
            path: path.to_owned(),
        })
    }

    pub(crate) fn create(path: &Path) -> Result<Self> {
        if path.exists() {
            return Err(Error::ExistingOutputDir {
                path: path.to_owned(),
            }
            .into());
        }
        fs::create_dir_all(path)?;

        Ok(Outdir {
            path: path.to_owned(),
        })
    }

    pub(crate) fn serialize_value<V: SerdeSerialize>(&self, name: &str, value: V) -> Result<()> {
        let pathname = format!("{}.mpk", name);
        let file = fs::File::create(self.join(pathname))?;
        value.serialize(&mut Serializer::new(file))?;
        Ok(())
    }

    pub(crate) fn deserialize_value<'a, V: SerdeDeserialize<'a>>(&self, name: &str) -> Result<V> {
        let pathname = format!("{}.mpk", name);
        Ok(V::deserialize(&mut Deserializer::new(fs::File::open(
            self.join(pathname),
        )?))?)
    }
}


