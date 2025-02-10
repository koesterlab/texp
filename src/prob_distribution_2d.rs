use bio::stats::LogProb;
use itertools::iproduct;
use kdtree::KdTree;
use num_traits::Float;
use ordered_float::OrderedFloat;
use serde_derive::{Deserialize, Serialize};
use std::collections::BTreeMap;

//--------------------ProbDistribution based on KdTree --------------------

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct EntryForKdtree {
    position: [f64; 2],
    prob: LogProb,
    trbl: [usize; 4], // top right bottom left
}

/// Datastructure for storing sample expression probability distributions. kdtree is a 2 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution2d {
    pub points: Vec<EntryForKdtree>,
    pub kdtree: KdTree<f64, usize, [f64; 2]>,
    max_prob_entry_position: usize,
    max_distance: f64,
    max_x: f64,
    max_y: f64,
    range_per_theta: BTreeMap<OrderedFloat<f64>, [f64; 2]>,
}

pub fn euclidean<T: Float>(a: &[T], b: &[T]) -> T {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| ((*x) - (*y)) * ((*x) - (*y)))
        .fold(T::zero(), ::std::ops::Add::add)
        .sqrt()
}

impl ProbDistribution2d {
    pub(crate) fn new() -> Self {
        ProbDistribution2d {
            points: Vec::<EntryForKdtree>::new(),
            kdtree: KdTree::new(2), // 2 dimensional kdTree
            max_prob_entry_position: usize::MAX,
            max_distance: 0.,
            max_x: 0.,
            max_y: 0.,
            range_per_theta: BTreeMap::new(),
        }
    }

    pub(crate) fn na() -> Self {
        ProbDistribution2d {
            points: Vec::<EntryForKdtree>::new(),
            kdtree: KdTree::new(2),
            max_prob_entry_position: usize::MAX,
            max_distance: 0.,
            max_x: 0.,
            max_y: 0.,
            range_per_theta: BTreeMap::new(),
        }
    }

    // pub(crate) fn len(&self) -> usize {
    //     self.points.len()
    // }

    pub(crate) fn get_max_prob(&self) -> LogProb {
        self.points[self.max_prob_entry_position].prob
    }

    // pub(crate) fn get_max_prob_keys(&self) -> [f64; 2] {
    //     self.points[self.max_prob_entry_position].position
    // }

    pub(crate) fn is_na(&self) -> bool {
        self.max_prob_entry_position == usize::MAX
    }

    pub(crate) fn insert_grid<F>(&mut self, mus: Vec<f64>, thetas: Vec<f64>, mut calc: F)
    where
        F: FnMut(f64, f64) -> LogProb,
    {
        self.max_y = thetas[thetas.len() - 1];
        self.max_x = mus[mus.len() - 1];
        // println!("max_x {:?}, max_y {:?}", self.max_x, self.max_y);
        let len_mus = mus.len();
        let len_thetas = thetas.len();
        for (j, i) in iproduct!(0..len_thetas, 0..len_mus) {
            let mu = mus[i];
            let theta = thetas[j];
            let own_pos = i + j * len_mus;
            let top = if j + 1 < len_thetas {
                own_pos + len_mus
            } else {
                usize::MAX
            };
            let right = if (i + 1) % len_mus != 0 {
                own_pos + 1
            } else {
                usize::MAX
            };
            let bottom = if j > 0 { own_pos - len_mus } else { usize::MAX };
            let left = if i % len_mus != 0 {
                own_pos - 1
            } else {
                usize::MAX
            };
            if mu.is_infinite() || theta.is_infinite() {
                println!("SHOULD NOT HAPPEN!!!");
                // println!("mu {:?}, theta {:?}, prob ---, top {:?}, right {:?}, bottom {:?}, left {:?}", mu, theta, top, right, bottom, left);
                let value: [f64; 2] = [mu, theta];
                let entry = EntryForKdtree {
                    position: value,
                    prob: LogProb::ln_zero(),
                    trbl: [top, right, bottom, left],
                };
                self.points.push(entry);
            } else {
                let prob = calc(mu, theta);
                // println!("mu {:?}, theta {:?}, prob {:?}, top {:?}, right {:?}, bottom {:?}, left {:?}", mu, theta, prob, top, right, bottom, left);
                self.insert(mu, theta, prob, top, right, bottom, left);
            }
        }
    }

    pub(crate) fn insert(
        &mut self,
        mu: f64,
        theta: f64,
        prob: LogProb,
        top: usize,
        right: usize,
        bottom: usize,
        left: usize,
    ) {
        let value: [f64; 2] = [mu, theta];
        let entry = EntryForKdtree {
            position: value,
            prob: prob,
            trbl: [top, right, bottom, left],
        };
        let entry_position = self.points.len();
        self.points.push(entry);

        if self.is_na() || self.get_max_prob() < prob {
            self.max_prob_entry_position = entry_position;
        }

        let distance = euclidean(&value, &[0.0, 0.0]);
        if distance.is_finite() {
            if distance > self.max_distance {
                self.max_distance = distance;
            }
            self.kdtree.add(value, entry_position).ok();
        }
        // if prob > LogProb::ln_zero(){
        //     if self.range_per_theta.contains_key(&OrderedFloat(theta)){
        //         let mu_range = self.range_per_theta.get(&OrderedFloat(theta)).unwrap();
        //             if mu < mu_range[0]{
        //                 self.range_per_theta.insert(OrderedFloat(theta), [mu, mu_range[1]]);
        //             } else if mu > mu_range[1]{
        //                 self.range_per_theta.insert(OrderedFloat(theta), [mu_range[0], mu]);
        //             }
        //     } else {
        //         self.range_per_theta.insert(OrderedFloat(theta), [mu, mu]);
        //     }

        // }
    }

    pub(crate) fn get(&self, value: &[f64]) -> LogProb {
        if self.is_na() {
            if value[0] == 0.0 {
                // mean 0
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            // for nearest conrer in neares_iter:
            // find square
            // if not in square: continue
            // println!("query mu {:?} theta {:?}",value[0], value[1]);
            let mut result = LogProb::ln_zero();
            if value[0] < 0. {
                // can happen in fold change calculation depending on constant c
                // println!("query negative, return 0");
                return result;
            }
            let mut nearest_iter = self.kdtree.iter_nearest(&value, &euclidean).unwrap();
            let nearest_corner_index = nearest_iter.next().unwrap().1;
            let nearest_point = &self.points[*nearest_corner_index];
            result = nearest_point.prob;

            if value[0] != nearest_point.position[0] {
                println!("query mu {:?} theta {:?}", value[0], value[1]);
                println!(
                    "nearest_point {:?} prob {:?}",
                    nearest_point, nearest_point.prob
                );
            }

            // let nearest_iter = self.kdtree.iter_nearest(&value, &euclidean).unwrap();
            // for nearest_corner in nearest_iter {
            //     let nearest_distance = nearest_corner.0;
            //     let nearest_corner_index = nearest_corner.1;
            //     let nearest_point = &self.points[*nearest_corner_index];
            // println!("nearest_point {:?} prob {:?}", nearest_point, nearest_point.prob);
            // let nearest_square = self.find_square_for_point(&value, *nearest_corner_index);
            // // println!("2d get neares_square {:?}", nearest_square);
            // if self.is_in_square(&value, &nearest_square){
            //     result = self.smoothing(&value, nearest_square);
            //     break;
            // }
            // }
            // println!("Result {:?}", result);
            return result;
        }
    }
}
