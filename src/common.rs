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
    prob_lower: LogProb,
    upper: N32,
    prob_upper: LogProb,
) -> LogProb {
    let len = upper - lower;
    (prob_lower + LogProb(f64::from(((upper - value) / len).ln())))
        .ln_add_exp(prob_upper + LogProb(f64::from(((value - lower) / len).ln())))
}

pub(crate) trait LikelihoodFunc {
    fn get(&self, values: &[N32]) -> LogProb;
    fn insert(&mut self, values: &[N32], prob: LogProb);
}

#[derive(Deserialize, Serialize, Debug)]
pub(crate) struct LikelihoodFunction<V> {
    points: BTreeMap<(N32, N32), LogProb>,
}

impl<V: LikelihoodFunc> LikelihoodFunc for LikelihoodFunction<V> {
    fn get(&self, values: &[N32]) -> LogProb {
        let upper = self.points.range(values[0]..).next();

        if let Some((upper, upper_prob)) = upper {
            if *upper == values[0] {
                upper_prob.get(&values[1..])
            } else {
                let lower = self.points.range(..=values[0]).last();
                if let Some((lower, lower_prob)) = lower {
                    interpolate_pmf(
                        values[0],
                        *lower,
                        lower_prob.get(&values[1..]),
                        *upper,
                        upper_prob.get(&values[1..]),
                    )
                } else {
                    LogProb::ln_zero()
                }
            }
        } else {
            LogProb::ln_zero()
        }
    }

    fn insert(&mut self, values: &[N32], prob: LogProb) {
        self.points.insert(values[0], 
    }
}

impl LikelihoodFunc for LikelihoodFunction<LogProb> {
    fn get(&self, values: &[N32]) -> LogProb {
        let upper = self.points.range(values[0]..).next();

        if let Some((upper, upper_prob)) = upper {
            if *upper == values[0] {
                *upper_prob
            } else {
                let lower = self.points.range(..=values[0]).last();
                if let Some((lower, lower_prob)) = lower {
                    interpolate_pmf(
                        values[0],
                        *lower,
                        *lower_prob,
                        *upper,
                        *upper_prob,
                    )
                } else {
                    LogProb::ln_zero()
                }
            }
        } else {
            LogProb::ln_zero()
        }
    }
}
