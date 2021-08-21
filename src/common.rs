use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use derefable::Derefable;
use itertools_num::linspace;
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
        linspace(mean / 5.0, mean, 10).rev().skip(1),
        linspace(mean, 5.0 * mean, 10),
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
    V: Ord + Eq,
{
    points: BTreeMap<V, LogProb>,
}

impl<V> ProbDistribution<V>
where
    V: Ord + Eq,
{
    pub(crate) fn get(&self, value: &V) -> LogProb {
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

    pub(crate) fn insert(&mut self, value: V, prob: LogProb) {
        self.points.insert(value, prob);
    }

    pub(crate) fn len(&self) -> usize {
        self.points.len()
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
        let file = fs::File::create(self.join(name).with_extension("mpk"))?;
        value.serialize(&mut Serializer::new(file))?;
        Ok(())
    }

    pub(crate) fn deserialize_value<'a, V: SerdeDeserialize<'a>>(&self, name: &str) -> Result<V> {
        Ok(V::deserialize(&mut Deserializer::new(fs::File::open(
            self.join(name).with_extension("mpk"),
        )?))?)
    }
}
