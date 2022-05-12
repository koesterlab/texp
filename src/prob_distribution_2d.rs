use anyhow::Result;
use bio::stats::LogProb;
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_derive::{Deserialize, Serialize};

//--------------------ProbDistribution based on KdTree --------------------

/// Datastructure for storing sample expression probability distributions. kdtree is a 2 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution2d {
    pub kdtree: KdTree<f64, LogProb, [f64; 2]>,
    max_prob_value: Option<[f64; 2]>,
    max_distance: f64,
    is_na: bool,
}

impl ProbDistribution2d {
    pub(crate) fn new() -> Self {
        ProbDistribution2d {
            kdtree: KdTree::new(2), // 2 dimensional kdTree
            max_prob_value: None,
            max_distance: 0.,
            is_na: true,
        }
    }

    pub(crate) fn na() -> Self {
        ProbDistribution2d {
            kdtree: KdTree::new(2),
            max_prob_value: None,
            max_distance: 0.,
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

        let distance = squared_euclidean(&value, &[0.0, 0.0]);
        if distance > self.max_distance {
            self.max_distance = distance;
        }        

        self.kdtree.add(value, prob)?;
        self.is_na = false;
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
            if squared_euclidean(&value, &[0.0, 0.0]) > self.max_distance * 10. {
                return LogProb::ln_zero();
            } else {
                return *self.kdtree.nearest(&value, 1, &squared_euclidean).unwrap()[0].1;
            }
        }
    }
}
