use std::collections::BTreeMap;

use bio::stats::LogProb;
use itertools_num::linspace;
use noisy_float::prelude::Float;
use noisy_float::types::N32;
use serde_derive::{Deserialize, Serialize};

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
    V: Ord + Eq
{
    points: BTreeMap<V, LogProb>,
}

impl<V> ProbDistribution<V> 
where
    V: Ord + Eq
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
