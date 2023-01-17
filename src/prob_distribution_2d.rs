use std::process::exit;
use std::collections::VecDeque;
use std::collections::HashMap;
use std::collections::BTreeMap;
use std::ops::Bound::Included;

use bio::stats::LogProb;
use itertools::iproduct;
use kdtree::KdTree;
use serde_derive::{Deserialize, Serialize};
use num_traits::Float;
use ordered_float::OrderedFloat;

use crate::common::{Square, Point, difference_to_big};

//--------------------ProbDistribution based on KdTree --------------------

const T: usize = 0;   // top
const R: usize = 1;   // right
const B: usize = 2;   // bottom
const L: usize = 3;   // left

#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct EntryForKdtree {
    position: [f64; 2],
    prob: LogProb,
    trbl: [usize; 4],   // top right bottom left
}

/// Datastructure for storing sample expression probability distributions. kdtree is a 2 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution2d {
    pub points: Vec<EntryForKdtree>,
    pub kdtree: KdTree<f64, usize, [f64; 2]>,
    max_prob_entry_position: usize,
    max_distance: f64,
    max_x:f64,
    max_y: f64,
    range_per_theta: BTreeMap<OrderedFloat<f64>, [f64; 2]>
}

pub fn euclidean<T: Float>(a: &[T], b: &[T]) -> T {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| ((*x) - (*y)) * ((*x) - (*y)))
        .fold(T::zero(), ::std::ops::Add::add).sqrt()
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

    pub(crate) fn len(&self) -> usize {
        self.points.len()
    }

    pub(crate) fn get_max_prob(&self) -> LogProb {
        self.points[self.max_prob_entry_position].prob
    }

    pub(crate) fn get_range_per_theta(&self, theta: f64) -> [f64; 2] {
        if self.range_per_theta.is_empty(){
            return [0.,10000.]; //TODO set variables as defaults
        }
        else if self.range_per_theta.contains_key(&OrderedFloat(theta)){
            return self.range_per_theta[&OrderedFloat(theta)]
        } else {
            // println!("theta {:?} not in map", theta);
            // let range = self.range_per_theta.range((Included(&OrderedFloat(theta/2.)), Included(&OrderedFloat(theta*2.))));
            // println!("range {:?}", range);
            // let mina = range.min_by(|x, y| ((x.0 - theta).abs()).cmp(&(y.0 - OrderedFloat(theta)).abs())).unwrap();
            // println!("mina {:?}", mina);
            *self.range_per_theta.range((Included(&OrderedFloat(theta/2.)), Included(&OrderedFloat(theta*2.))))
                .min_by(|x, y|((x.0 - theta).abs()).cmp(&(y.0 - OrderedFloat(theta)).abs()))
                .unwrap().1
        }
    }


    pub(crate) fn is_na(&self) -> bool {
        self.max_prob_entry_position == usize::MAX
    }

    pub(crate) fn check_if_different(p1: &Point, p2: &[f64; 2], reference: &[f64; 2]) -> bool {
        if (p1.x - p2[0]).abs() / reference[0] > 1e-6 {
            return true;
        }
        if (p1.y - p2[1]).abs() / reference[1] > 1e-6 {
            return true;
        }
        return false;
    }

    pub(crate) fn check_if_same(a:f64, b: f64) -> bool {
        if a == 0. && b == 0. {
            return true;
        }
        // println!("a {:?},b {:?}, diff {:?}, max {:?},result {:?}", a, b,(a - b).abs(), a.abs().max(b.abs()),  (a - b).abs() / a.abs().max(b.abs()) );
        if (a - b).abs() / a.abs().max(b.abs()) < 1e-6 {
            return true;
        }
        return false;
    }

    pub(crate) fn insert_grid<F>(&mut self, mut mus: Vec<f64>, mut thetas: Vec<f64>, calc: F) where F: Fn(f64, f64) -> LogProb {
        self.max_y = thetas[thetas.len()-1];
        self.max_x = mus[mus.len()-1];
        let len_mus = mus.len();
        let len_thetas = thetas.len();
        for (j, i) in iproduct!(0..len_thetas, 0..len_mus) {
            let mu = mus[i];
            let theta = thetas[j];
            let own_pos = i + j * len_mus;
            let top = if j + 1 < len_thetas {own_pos + len_mus} else {usize::MAX};
            let right = if (i + 1) % len_mus != 0 {own_pos + 1} else {usize::MAX};
            let bottom = if j > 0 {own_pos - len_mus} else {usize::MAX};
            let left = if i % len_mus != 0 {own_pos - 1} else {usize::MAX};
            if mu.is_infinite() || theta.is_infinite() {
                println!("SHOULD NOT HAPPEN!!!");
                // println!("mu {:?}, theta {:?}, prob ---, top {:?}, right {:?}, bottom {:?}, left {:?}", mu, theta, top, right, bottom, left);
                let value: [f64; 2] = [mu, theta];
                let entry = EntryForKdtree{position: value, prob: LogProb::ln_zero(), trbl: [top, right, bottom, left]};
                self.points.push(entry);
            } else {
                let prob = calc(mu, theta);
                // println!("mu {:?}, theta {:?}, prob {:?}, top {:?}, right {:?}, bottom {:?}, left {:?}", mu, theta, prob, top, right, bottom, left);
                self.insert(mu, theta, prob, top, right, bottom, left);
            }
        }
        
        let mut queue = VecDeque::<Square>::new();
        for (j, i) in iproduct!(0..len_thetas - 1, 0..len_mus - 1) {
            queue.push_back( Square { top_left: i + j * len_mus + len_mus,
                top_right: i + j * len_mus + len_mus + 1,
                bot_left: i + j * len_mus,
                bot_right: i + j * len_mus + 1,
            });
        }

        while queue.len() > 0 {
            let mut debugPrint = false;
            let square = queue.pop_front().unwrap();
            let new_point_candidates = self.new_points(&square);
            let mut calced_values = Vec::<Option<LogProb>>::new();
            let mut estimated_values = Vec::<Option<LogProb>>::new();
            let mut new_point_indices = Vec::<usize>::new();
            let mut counter = self.len();
            for p in &new_point_candidates {
                if p.x.is_infinite() || p.y.is_infinite(){
                    calced_values.push(Some(LogProb::ln_zero()));
                    estimated_values.push(Some(LogProb::ln_zero()));
                    new_point_indices.push(counter);
                    counter += 1;
                    continue;
                }
                let nearest = self.kdtree.nearest(&[p.x, p.y], 1, &euclidean).unwrap()[0];
                let distance_to_nearest = nearest.0;
                // println!("counter {:?}, p.x {:?}, p.y {:?}, dist {:?}, nearst {:?}", counter, p.x, p.y, distance_to_nearest, self.points[*nearest.1]);
                if ProbDistribution2d::check_if_different(p, &self.points[*nearest.1].position, &self.points[square.bot_left].position) {
                        let prob = calc(
                        p.x,  // mu_ik
                        p.y,  // *theta_i
                    );
                    // println!("prob {:?}, p.x {:?}, p.y {:?}", prob, p.x, p.y);
                    calced_values.push(Some(prob));
                    estimated_values.push(Some(self.get(&[p.x, p.y])));
                    new_point_indices.push(counter);
                    counter += 1;
                } else {
                    calced_values.push(None);
                    estimated_values.push(None);
                    new_point_indices.push(*nearest.1);
                }
            }
            // println!("calced {:?}, est {:?}, new_p_index {:?}", calced_values, estimated_values, new_point_indices);
            let mut need_iteration = false;
            for (estimated_value, calculated_value) in estimated_values.iter().zip(calced_values.iter()) {
                if estimated_value.is_none() {
                    continue;
                }
                let prob_dist_max = self.get_max_prob();
                if difference_to_big(estimated_value.unwrap(), calculated_value.unwrap(), prob_dist_max) {
                    need_iteration = true;
                    break;
                }
            }
            if !calced_values[0].is_none() {
                let point = new_point_candidates[0];
                self.insert(point.x, point.y, calced_values[0].unwrap(), new_point_indices[2], square.bot_right, usize::MAX, square.bot_left);
                self.points[square.bot_right].trbl[L] = new_point_indices[0];
            self.points[square.bot_left].trbl[R] = new_point_indices[0];
            if calced_values[2].is_none() && new_point_indices[2] != usize::MAX {
                self.points[new_point_indices[2]].trbl[B] = new_point_indices[0];
            }
                // auch für die anderen Punkte testen, ob Sie schon existieren und dann top/left/... von denen entsprechend setzen
            }   
            
            if !calced_values[1].is_none() {
                let point = new_point_candidates[1];
                self.insert(point.x, point.y, calced_values[1].unwrap(), square.top_left, new_point_indices[2], square.bot_left, usize::MAX);
                self.points[square.top_left].trbl[B] = new_point_indices[1];
                self.points[square.bot_left].trbl[T] = new_point_indices[1];
                if calced_values[2].is_none() && new_point_indices[2] != usize::MAX {
                    self.points[new_point_indices[2]].trbl[L] = new_point_indices[1];
                }
            }
            
            if !calced_values[2].is_none() {
                let point = new_point_candidates[2];
                self.insert(point.x, point.y, calced_values[2].unwrap(), new_point_indices[4], new_point_indices[3], new_point_indices[0], new_point_indices[1]);
                if new_point_indices[0] != usize::MAX {
                    self.points[new_point_indices[0]].trbl[T] = new_point_indices[2];
                }
                if new_point_indices[1] != usize::MAX {
                    self.points[new_point_indices[1]].trbl[R] = new_point_indices[2];
                }
                if new_point_indices[3] < self.len() {
                    self.points[new_point_indices[3]].trbl[L] = new_point_indices[2];
                }
                if new_point_indices[4] < self.len() {
                    self.points[new_point_indices[4]].trbl[B] = new_point_indices[2];
                }
            }

            if !calced_values[3].is_none() {
                let point = new_point_candidates[3];
                self.insert(point.x, point.y, calced_values[3].unwrap(), square.top_right, usize::MAX, square.bot_right, new_point_indices[2]);
                self.points[square.top_right].trbl[B] = new_point_indices[3];
                self.points[square.bot_right].trbl[T] = new_point_indices[3];
                if calced_values[2].is_none() && new_point_indices[2] != usize::MAX {
                    self.points[new_point_indices[2]].trbl[R] = new_point_indices[3];
                }
            }
            
            if !calced_values[4].is_none() {
                let point = new_point_candidates[4];
                self.insert(point.x, point.y, calced_values[4].unwrap(), usize::MAX, square.top_right, new_point_indices[2], square.top_left);
                self.points[square.top_right].trbl[L] = new_point_indices[4];
                self.points[square.top_left].trbl[R] = new_point_indices[4];
                if calced_values[2].is_none() && new_point_indices[2] != usize::MAX {
                    self.points[new_point_indices[2]].trbl[T] = new_point_indices[4];
                }
            }
            
            if self.len() > 10000 {
                need_iteration = false;
            }
            if need_iteration {
                queue.push_back(Square { top_left: square.top_left,
                    top_right: new_point_indices[4],
                    bot_left: new_point_indices[1],
                    bot_right: new_point_indices[2],
                });
                queue.push_back(Square { top_left: new_point_indices[4],
                    top_right: square.top_right,
                    bot_left: new_point_indices[2],
                    bot_right: new_point_indices[3],
                });
                queue.push_back(Square { top_left: new_point_indices[1],
                    top_right: new_point_indices[2],
                    bot_left: square.bot_left,
                    bot_right: new_point_indices[0],
                });
                queue.push_back(Square { top_left: new_point_indices[2],
                    top_right: new_point_indices[3],
                    bot_left: new_point_indices[0],
                    bot_right: square.bot_right,
                });
            }
        }
    }


    pub(crate) fn new_points(&mut self, square: &Square) -> Vec<Point> {
        let mut result = Vec::<Point>::new();
    
        let new_point = |p1_pos: usize, p2_pos: usize| {
            let p1 = self.points[p1_pos].position;
            let p2 = self.points[p2_pos].position;
            let new_coord = |left: f64, right: f64| {
                let middle; //= 0.;
                if right.is_finite() {
                    middle = left / 2. + right / 2.;
                } else {
                    middle = 10. * left;
                }
                return middle;
            };
            return Point{x:new_coord(p1[0], p2[0]), y:new_coord(p1[1], p2[1])};
        };
    
        result.push(new_point(square.bot_left, square.bot_right));
        result.push(new_point(square.bot_left, square.top_left));
        result.push(new_point(square.bot_left, square.top_right));
        result.push(new_point(square.bot_right, square.top_right));
        result.push(new_point(square.top_left, square.top_right));
    
        return result;
    }

    pub(crate) fn insert(&mut self, mu: f64, theta: f64, prob: LogProb, top: usize, right: usize, bottom: usize, left: usize) {
        // if mu > 47. && mu < 51. {
        //     println!("insert: mu {:?}, theta {:?}, prob {:?}", mu, theta, prob.exp());
        // }
        let value: [f64; 2] = [mu, theta];
        let entry = EntryForKdtree{position: value, prob: prob, trbl: [top, right, bottom, left]};
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
        if prob > LogProb::ln_zero(){
            if self.range_per_theta.contains_key(&OrderedFloat(theta)){
                let mu_range = self.range_per_theta.get(&OrderedFloat(theta)).unwrap();
                    if mu < mu_range[0]{
                        self.range_per_theta.insert(OrderedFloat(theta), [mu, mu_range[1]]);
                    } else if mu > mu_range[1]{
                        self.range_per_theta.insert(OrderedFloat(theta), [mu_range[0], mu]);
                    }
            } else {
                self.range_per_theta.insert(OrderedFloat(theta), [mu, mu]);
            }

        }
    }

    pub(crate) fn get(&self, value: &[f64]) -> LogProb {
        if self.is_na() {
            if value[0] == 0.0 {
                // mean 0
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else { // for nearest conrer in neares_iter:
            // find square
            // if not in square: continue
                let mut result = LogProb::ln_zero();
                let nearest_iter = self.kdtree.iter_nearest(&value, &euclidean).unwrap();
                for nearest_corner in nearest_iter {
                    let nearest_distance = nearest_corner.0;
                    let nearest_corner_index = nearest_corner.1;
                    let nearest_square = self.find_square_for_point(&value, *nearest_corner_index);
                    // println!("2d get neares_square {:?}", nearest_square);
                    if self.is_in_square(&value, &nearest_square){
                        result = self.smoothing(&value, nearest_square);
                        break;
                    }                   
                }
            return result;
        }
    }

    pub(crate) fn push_if(vector: &mut Vec<usize>, value: usize, position: &[f64; 2]) {
        if vector.contains(&value) {
            println!("Endless loop detected");
            exit(1);
        }
        if position[0].is_infinite() || position[1].is_infinite() {
            return;
        }
        vector.push(value);
    }

    pub(crate) fn find_square_for_point(&self, point: &[f64], nearest_point_index: usize) -> Vec<usize> {
        // Starte beim nearest_point_index Punkt, bestimme eine Richtung abhängig vom point -> TRBL und Richtungsänderung CW +1 oder CCW -1
        // Gehe eine Richtung bis abzweigbar ist, dann Richtung + Richtungsänderung und wiederholen bis man am Ursprung oder am Rand ist
        // println!("-----Start find_square_for_point {:?}", point);
        let nearest_point = &self.points[nearest_point_index];
        let mut result = Vec::<usize>::new();
        let mut direction;
        let mut direction_change;
        if point[0] < nearest_point.position[0] {
            direction = L;
        } else {
            direction = R;
        }
        if point[1] < nearest_point.position[1] {
            direction_change = if direction == L {3} else {1};  // 3 == -1 mod 4
        } else if ProbDistribution2d::check_if_same(point[1], nearest_point.position[1]) 
                && ProbDistribution2d::check_if_same(point[1], self.max_y) { //point[1] == nearest_point.position[1] && point[1] == self.max_y {
            direction_change = if direction == L {3} else {1};
        } else {
            direction_change = if direction == L {1} else {3};  // 3 == -1 mod 4
        }

        if point[0] >= self.max_x {
            // println!("X out of bound");
            ProbDistribution2d::push_if(&mut result, nearest_point_index, &nearest_point.position);
            let mut current_point_index = nearest_point.trbl[0];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                current_point_index = self.points[current_point_index].trbl[0];
                // println!("curr {:?}", current_point_index);
                if current_point_index != usize::MAX {
                    ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                }
            }
            current_point_index = nearest_point.trbl[3];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
            }
            current_point_index = nearest_point.trbl[2];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                current_point_index = self.points[current_point_index].trbl[2];
                // println!("curr {:?}", current_point_index);
                if current_point_index != usize::MAX {
                    ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                }
            }
            return result;
        }

        if point[1] >= self.max_y {
            // println!("Y out of bound");
            ProbDistribution2d::push_if(&mut result, nearest_point_index, &nearest_point.position);
            let mut current_point_index = nearest_point.trbl[1];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                current_point_index = self.points[current_point_index].trbl[1];
                // println!("curr {:?}", current_point_index);
                if current_point_index != usize::MAX {
                    ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                }
            }
            current_point_index = nearest_point.trbl[2];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
            }
            current_point_index = nearest_point.trbl[3];
            // println!("curr {:?}", current_point_index);
            if current_point_index != usize::MAX {
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                current_point_index = self.points[current_point_index].trbl[3];
                // println!("curr {:?}", current_point_index);
                if current_point_index != usize::MAX {
                    ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
                }
            }
            return result;
        }

        // println!("neare {:?}, near_p {:?}, dir {:?}, dir_change {:?}", nearest_point_index, nearest_point, direction, direction_change);
        ProbDistribution2d::push_if(&mut result, nearest_point_index, &nearest_point.position);
        let mut current_point_index = nearest_point.trbl[direction];
        if current_point_index == usize::MAX {
            direction = (direction + direction_change) % 4;
            direction_change = (direction_change + 2) % 4;
            current_point_index = nearest_point.trbl[direction];
            // println!("CHANGE dir {:?}, dir_change {:?}, curr_in {:?}", direction, direction_change, current_point_index);
        }
        while current_point_index != nearest_point_index {
            ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
            let next_direction = (direction + direction_change) % 4;
            // println!("curr {:?}, point {:?}, dir {:?}, dir_change {:?}, next_dir {:?}, value {:?}", current_point_index, self.points[current_point_index], direction, direction_change, next_direction, self.points[current_point_index].trbl[next_direction]);
            while self.points[current_point_index].trbl[next_direction] == usize::MAX {
                current_point_index = self.points[current_point_index].trbl[direction];
                // println!("curr {:?}, point {:?}, dir {:?}, dir_change {:?}, next_dir {:?}, value {:?}", current_point_index, self.points[current_point_index], direction, direction_change, next_direction, self.points[current_point_index].trbl[next_direction]);
                if current_point_index == nearest_point_index {
                    break;
                }
                ProbDistribution2d::push_if(&mut result, current_point_index, &self.points[current_point_index].position);
            }
            if current_point_index == nearest_point_index {
                break;
            }
            direction = next_direction;
            current_point_index = self.points[current_point_index].trbl[direction];
        }

        return result;
    }


    pub(crate) fn is_in_square(&self, value: &[f64], square: &Vec<usize> ) -> bool {
        let mut min_x = f64::MAX;
        let mut min_y = f64::MAX;
        let mut max_x = f64::MIN;
        let mut max_y = f64::MIN;
        for index in square.iter(){
            let p = &self.points[*index];
            min_x = min_x.min(p.position[0]);
            min_y = min_y.min(p.position[1]);
            max_x = max_x.max(p.position[0]);
            max_y = max_y.max(p.position[1]);
        }
        if value[0] >= min_x && value[0] <= max_x && value[1] >= min_y && value[1] <= max_y {
            return true;
        } else {
            return false;
        }
    }


    pub(crate) fn smoothing(&self, value: &[f64], point_indices: Vec<usize>) -> LogProb {
        let c = 1.;
        let mut min_x = f64::MAX;
        let mut min_y = f64::MAX;
        let mut max_x = f64::MIN;
        let mut max_y = f64::MIN;
        let filtering = |x: &usize| {
            let p = &self.points[*x];
            let x1 = p.position[0] - value[0];
            let y1 = p.position[1] - value[1];
            let denominator_left = x1 * x1 + y1 * y1;
            if denominator_left.is_infinite() {
                return None;
            } else {
                min_x = min_x.min(p.position[0]);
                min_y = min_y.min(p.position[1]);
                max_x = max_x.max(p.position[0]);
                max_y = max_y.max(p.position[1]);
                return Some(p);
            }
        };
        
        let nearest_vec = point_indices.iter().filter_map(filtering).collect::<Vec<_>>();
        // println!("nearest_vec {:?}", nearest_vec);
        if nearest_vec.len() == 0 {
            return LogProb::ln_zero();
        }
        let mut result = LogProb::ln_zero();
        let mut scaling_sum = 0.;
        let mut distances_x = vec![0.; nearest_vec.len()];
        let mut distances_y = vec![0.; nearest_vec.len()];

        for i in 0..nearest_vec.len() {
            distances_x[i] = (value[0] - nearest_vec[i].position[0]).abs();//
            distances_y[i] = (value[1] - nearest_vec[i].position[1]).abs() ;//
        }

        for i in 0..nearest_vec.len() {
            let mut edge_scale_x = (1. - distances_x[i] / (max_x - min_x)).powi(8);
            let mut edge_scale_y = (1. - distances_y[i] / (max_y - min_y)).powi(8);
            edge_scale_x = (-1. / edge_scale_x).exp();
            edge_scale_y = (-1. / edge_scale_y).exp();

            scaling_sum += 1. * edge_scale_x * edge_scale_y; 
            result = result.ln_add_exp(nearest_vec[i].prob    + LogProb(edge_scale_x.ln()) + LogProb(edge_scale_y.ln()));
        } 
        if result == LogProb::ln_zero(){
            return result;
        }

        result = result - LogProb(scaling_sum.ln());
        let distance_of_value =  euclidean(&value, &[0.0, 0.0]);
        if self.max_distance < distance_of_value { //is_outlier
            result = result - LogProb((distance_of_value - self.max_distance).ln());  
            println!("Is outlier. New result {:?}", result.exp())      
        }
        return result;
    }
}
