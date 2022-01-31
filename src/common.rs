use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use derefable::Derefable;
use derive_new::new;
use itertools_num::linspace;
use noisy_float::NoisyFloat;
use noisy_float::prelude::Float;
use noisy_float::types::N32;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};

use crate::errors::Error;

pub(crate) fn window(mean: f64) -> (impl Iterator<Item = f64>, impl Iterator<Item = f64>) {
    // TODO: think about larger steps, binary search etc. to optimize instead of just having a fixed number of steps.
    (
        linspace(mean / 5.0, mean, 15).rev().skip(1),
        linspace(mean, 5.0 * mean, 25),
    )
}

pub(crate) fn interpolate_pmf(
    value: N32,
    lower: N32,
    upper: N32,
    prob_lower: LogProb,
    prob_upper: LogProb,
) -> LogProb {
    let len = upper - lower;
    (prob_lower + LogProb(f64::from(((upper - value) / len).ln())))
        .ln_add_exp(prob_upper + LogProb(f64::from(((value - lower) / len).ln())))
}

#[derive(Deserialize, Serialize, Debug, Default)]
pub(crate) struct ProbDistribution<V>
where
    V: Ord + Eq + Copy + DistributionValue,
{
    pub points: BTreeMap<V, LogProb>,
    max_prob_value: Option<V>,
    is_na: bool,
}

impl ProbDistribution<Mean> {
    pub(crate) fn normalize(&mut self) {
        let marginal = LogProb::ln_trapezoidal_integrate_grid_exp(
    |i, value| {
                *self.points.get(&Mean::new(value)).unwrap()
            }, 
        &self.points.keys().map(|value| **value).collect::<Vec<_>>()
        );

        for prob in self.points.values_mut() {
            *prob = *prob - marginal
        }
    }
}

impl<V> ProbDistribution<V>
where
    V: Ord + Eq + Copy + DistributionValue,
{
    pub(crate) fn na() -> Self {
        ProbDistribution {
            points: BTreeMap::default(),
            max_prob_value: None,
            is_na: true,
        }
    }

    pub(crate) fn get(&self, value: &V) -> LogProb {
        if self.is_na {
            if value.is_zero() {
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            let upper = self.points.range(value..).next();

            if let Some((upper, upper_prob)) = upper {
                if upper == value {
                    *upper_prob
                } else {
                    let lower = self.points.range(..=value).last();
                    if let Some((lower, lower_prob)) = lower {
                        // TODO interpolate
                        LogProb(*lower_prob.ln_add_exp(*upper_prob) - 2.0_f64.ln())
                    } else {
                        LogProb::ln_zero()
                    }
                }
            } else {
                LogProb::ln_zero()
            }
        }
    }

    pub(crate) fn insert(&mut self, value: V, prob: LogProb) {
        self.points.insert(value, prob);

        if self.max_prob_value.is_none() || self.points.get(&self.max_prob_value.unwrap()).unwrap() < &prob {
            self.max_prob_value = Some(value);
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.points.len()
    }

    pub(crate) fn get_max_prob_value(&self) -> Option<V> {
        self.max_prob_value
    }
}

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

pub(crate) trait DistributionValue {
    fn is_zero(&self) -> bool;
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy, Deserialize, Serialize, new, Default)]
pub(crate) struct MeanDispersionPair {
    mean: N32,
    Dispersion: N32,
}


impl DistributionValue for MeanDispersionPair {
    fn is_zero(&self) -> bool {
        self.mean == N32::new(0.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy, Deserialize, Serialize, new, Default, Derefable)]
pub(crate) struct Mean(#[deref] N32);


impl DistributionValue for Mean {
    fn is_zero(&self) -> bool {
        **self == N32::new(0.0)
    }
}


#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy, Deserialize, Serialize, new, Default, Derefable)]
pub(crate) struct Log2FoldChange(#[deref] N32);


impl DistributionValue for Log2FoldChange {
    fn is_zero(&self) -> bool {
        **self == N32::new(0.0)
    }
}