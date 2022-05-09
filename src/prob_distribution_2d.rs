use std::collections::BTreeMap;

use anyhow::Result;
use bio::stats::LogProb;
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use noisy_float::prelude::Float;
use noisy_float::types::N32;
use noisy_float::types::N64;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};

//--------------------ProbDistribution based on KdTree --------------------

/// Datastructure for storing sample expression probability distributions. kdtree is a 2 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution2d {
    points: BTreeMap<[N64; 2], LogProb>,
    pub kdtree: KdTree<f64, LogProb, [f64; 2]>,
    max_prob_value: Option<[f64; 2]>,
    is_na: bool,
}

impl ProbDistribution2d {
    pub(crate) fn new() -> Self {
        ProbDistribution2d {
            points: BTreeMap::default(),
            kdtree: KdTree::new(2), // 2 dimensional kdTree
            max_prob_value: None,
            is_na: true,
        }
    }

    pub(crate) fn na() -> Self {
        ProbDistribution2d {
            points: BTreeMap::default(),
            kdtree: KdTree::new(2),
            max_prob_value: None,
            is_na: true,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.kdtree.size()
    }

    pub(crate) fn get_max_prob_value(&self) -> [f64; 2] {
        self.max_prob_value.unwrap()
    }

    pub(crate) fn insert(&mut self, mu: f64, theta: f64, prob: LogProb) -> Result<()> {
        let value: [f64; 2] = [mu, theta];
        if self.is_na || self.kdtree.nearest(&value, 1, &squared_euclidean).unwrap()[0].1 < &prob {
            self.max_prob_value = Some(value);
        }

        self.kdtree.add(value, prob)?;
        self.is_na = false;

        let value: [N64; 2] = [N64::new(mu), N64::new(theta)];
        self.points.insert(value, prob);
        Ok(())
    }

    pub(crate) fn get(&self, value: &[f64]) -> LogProb {
        if self.is_na {
            if value[0] == 0.0 {
                // mean 0
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            *self.kdtree.nearest(&value, 1, &squared_euclidean).unwrap()[0].1
        }
    }
}
