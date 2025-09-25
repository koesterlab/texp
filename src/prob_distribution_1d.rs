use std::cmp::min;
use std::collections::BTreeMap;

// use anyhow::Result;
use bio::stats::LogProb;
use noisy_float::types::N64;
use serde_derive::{Deserialize, Serialize};

/// Datastructure for storing group expression probability distributions and fold change distributions. kdtree is a 1 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution1d {
    // pub points: BTreeMap<N64, (LogProb, f64, LogProb)>,
    pub points: BTreeMap<N64, LogProb>,
    max_prob_position: Option<f64>,
    max_prob: Option<LogProb>,
    is_na: bool,
}

impl ProbDistribution1d {
    pub(crate) fn new() -> Self {
        ProbDistribution1d {
            points: BTreeMap::default(),
            max_prob_position: None,
            max_prob: None,
            is_na: true,
        }
    }

    #[allow(unused)]
    pub(crate) fn na() -> Self {
        ProbDistribution1d {
            points: BTreeMap::default(),
            max_prob_position: None,
            max_prob: None,
            is_na: true,
        }
    }

    #[allow(unused)]
    pub(crate) fn len(&self) -> usize {
        self.points.len()
    }

    #[allow(unused)]
    pub(crate) fn get_max_prob_position(&self) -> f64 {
        self.max_prob_position.unwrap()
    }

    pub(crate) fn get_max_prob(&self) -> LogProb {
        self.max_prob.unwrap()
    }

    pub(crate) fn insert(&mut self, value: f64, prob: LogProb) {
        // println!("value {:?}, prob {:?}, size {:?}", value, prob, self.points.len());
        if value == f64::INFINITY || prob == LogProb::from(f64::INFINITY) {
            println!("value inf, prob {:?}, size {:?}", prob, self.points.len());
        }
        let value2 = [value];
        if self.is_na || self.max_prob.unwrap() < prob {
            self.max_prob_position = Some(value2[0]);
            self.max_prob = Some(prob);
        }
        self.is_na = false;
        self.points.insert(N64::new(value), prob);
    }

    pub(crate) fn get(&self, value: f64) -> LogProb {
        if self.is_na || value == f64::INFINITY {
            if value == 0.0 {
                // mean or fold change 0
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            //  println!("query f {:?} ",value);
            if let Some(prob) = self.points.get(&N64::new(value)) {
                return *prob;
            } else {
                return LogProb::ln_zero();
            }
        }
    }

    // pub(crate) fn normalize(&mut self) -> LogProb {
    //     if self.is_na {
    //         return LogProb::ln_one();
    //     }
    //     let density = |_, value| self.get(value);
    //     let marginals = self
    //         .points
    //         .keys()
    //         .map(|value| *value)
    //         .collect::<Vec<_>>()
    //         .windows(2)
    //         .map(|x| LogProb::ln_simpsons_integrate_exp(density, x[0].raw(), x[1].raw(), 3))
    //         .collect::<Vec<_>>();
    //     let marginal = LogProb::ln_sum_exp(&marginals);
    //     if marginal != LogProb::ln_zero() {
    //         for (prob, _, _) in self.points.values_mut() {  // d1, d2
    //             *prob = *prob - marginal //Logspace / -> -
    //         }
    //         self.max_prob = Some(self.max_prob.unwrap() - marginal);
    //     }
    //     marginal
    // }
}
