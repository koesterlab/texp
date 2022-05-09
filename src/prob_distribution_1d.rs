use std::cmp::min;
use std::collections::BTreeMap;

use anyhow::Result;
use bio::stats::LogProb;
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use noisy_float::types::N64;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use serde_derive::{Deserialize, Serialize};

/// Datastructure for storing group expression probability distributions and fold change distributions. kdtree is a 1 dimensional kdTree with data = probability in LogProb.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct ProbDistribution1d {
    points: BTreeMap<N64, (LogProb, f64, LogProb)>,
    // pub kdtree: KdTree<f64, LogProb, [f64; 1]>,
    max_prob_value: Option<f64>,
    is_na: bool,
}

impl ProbDistribution1d {
    pub(crate) fn new() -> Self {
        ProbDistribution1d {
            points: BTreeMap::default(),
            // kdtree: KdTree::with_capacity(1, 32), // 2 dimensional kdTree
            max_prob_value: None,
            is_na: true,
        }
    }

    pub(crate) fn na() -> Self {
        ProbDistribution1d {
            points: BTreeMap::default(),
            // kdtree: KdTree::new(1),
            max_prob_value: None,
            is_na: true,
        }
    }

    pub(crate) fn len(&self) -> usize {
        // self.kdtree.size()
        self.points.len()
    }

    pub(crate) fn get_max_prob_value(&self) -> f64 {
        self.max_prob_value.unwrap()
    }

    fn calc_directions(
        x1: f64,
        y1: LogProb,
        s1: f64,
        x2: f64,
        y2: LogProb,
        s2: f64,
    ) -> (f64, LogProb, f64) {
        let y1_ = f64::from(y1).exp();
        let y2_ = f64::from(y2).exp();
        let d1 = 0.; // (x1/(x1*x1+ y1_*y1_).sqrt() + x2/(x2*x2+ y2_*y2_).sqrt()); // is f64

        let x1_ = LogProb(x1.abs().ln());
        let x2_ = LogProb(x2.abs().ln());
        let mut left = LogProb::ln_zero();
        if !(y1 == LogProb::ln_zero() && x1 == 0.) {
            left = y1 - LogProb(f64::from((x1_ + x1_).ln_add_exp(y1 + y1)) / 2.);
        }
        let mut right = LogProb::ln_zero();
        if !(y2 == LogProb::ln_zero() && x2 == 0.) {
            right = y2 - LogProb(f64::from((x2_ + x2_).ln_add_exp(y2 + y2)) / 2.);
        }
        let mut d2 = LogProb::ln_zero();
        let mut new_s = s1;
        if (s1 == s2) {
            d2 = left.ln_add_exp(right);
        } else if (left < right) {
            new_s = s2;
            d2 = right.ln_sub_exp(left);
        } else {
            d2 = left.ln_sub_exp(right);
        }
        (d1, d2, new_s)
    }

    fn calc_sign_params(
        lower_prob: LogProb,
        prob: LogProb,
        upper_prob: LogProb,
    ) -> (f64, LogProb, f64, LogProb) {
        let mut sign1 = 1.;
        let mut param1 = LogProb::ln_zero();
        // Upper and lower bound is there
        if lower_prob < prob {
            sign1 = -1.;
            param1 = prob.ln_sub_exp(lower_prob)
        } else {
            param1 = lower_prob.ln_sub_exp(prob)
        }
        let mut sign2 = 1.;
        let mut param2 = LogProb::ln_zero();
        if upper_prob < prob {
            // println!("panic 2");
            sign2 = -1.;
            param2 = prob.ln_sub_exp(upper_prob)
        } else {
            param2 = upper_prob.ln_sub_exp(prob)
        }
        (sign1, param1, sign2, param2)
    }

    pub(crate) fn insert(&mut self, value: f64, prob: LogProb) -> Result<()> {
        if value == f64::INFINITY {
            println!("value inf, prob {:?}, size {:?}", prob, self.points.len());
        }
        let value2 = [value];
        // if self.is_na || self.kdtree.nearest(&value2, 1, &squared_euclidean).unwrap()[0].1 < &prob {
        //     self.max_prob_value = Some(value2[0]);
        // }
        // self.kdtree.add(value2, prob)?;

        if self.is_na
            || self
                .points
                .get(&N64::new(self.max_prob_value.unwrap()))
                .unwrap()
                .0
                < prob
        {
            self.max_prob_value = Some(value2[0]);
        }
        self.is_na = false;
        let mut d1 = 1.;
        let mut d2 = LogProb::ln_zero();
        let mut un_d1 = 1.;
        let mut un_d2 = LogProb::ln_zero();
        let mut lp_d1 = 1.;
        let mut lp_d2 = LogProb::ln_zero();

        {
            // let upper = self.points.range(N64::new(value)..).next();
            // let lower = self.points.range(..=N64::new(value)).last();
            let mut upper_it = self.points.range(N64::new(value)..);
            let upper = upper_it.next();
            let upper_next = upper_it.next();
            let mut lower_it = self.points.range(..=N64::new(value)).rev();
            let lower = lower_it.next();
            let lower_prev = lower_it.next();

            if let Some((upper, (upper_prob, ud1, ud2))) = upper {
                if let Some((lower, (lower_prob, ld1, ld2))) = lower {
                    let (sign1, param1, sign2, param2) =
                        ProbDistribution1d::calc_sign_params(*lower_prob, prob, *upper_prob);
                    let (local_d1, local_d2, local_s) = ProbDistribution1d::calc_directions(
                        lower.raw() - value,
                        param1,
                        sign1,
                        upper.raw() - value,
                        param2,
                        sign2,
                    );
                    // d1 = local_d1;
                    // d1 was changed to s1 Vorzeichen der LogProb, um negative LogProbs speichern zu können
                    d1 = local_s;
                    d2 = local_d2;

                    if let Some((upper_next, (upper_next_prob, u_next_d1, u_next_d2))) = upper_next
                    {
                        un_d1 = *u_next_d1;
                        un_d2 = *u_next_d2;
                    }
                    if let Some((lower_prev, (lower_prev_prob, l_prev_d1, l_prev_d2))) = lower_prev
                    {
                        lp_d1 = *l_prev_d1;
                        lp_d2 = *l_prev_d2;
                    }
                }
            }
            //     If only one bound or neither upper nor lower bound  are there
            //     d1 and d2 default to 0
        }

        {
            let upper = self.points.range_mut(N64::new(value)..).next();
            if let Some((upper, (upper_prob, ud1, ud2))) = upper {
                let (sign1, param1, sign2, param2) =
                    ProbDistribution1d::calc_sign_params(prob, *upper_prob, un_d2);
                let (local_d1, local_d2, local_s) = ProbDistribution1d::calc_directions(
                    // lower.raw()-value, lower_prob.ln_sub_exp(prob), upper.raw()-value, upper_prob.ln_sub_exp(prob));
                    // value-upper.raw(), prob.ln_sub_exp(*upper_prob), un_d1-upper.raw(), un_d2.ln_sub_exp(*upper_prob));
                    value - upper.raw(),
                    param1,
                    sign1,
                    un_d1 - upper.raw(),
                    param2,
                    sign2,
                );
                *ud1 = local_s;
                *ud2 = local_d2;
            }
        }

        {
            let lower = self.points.range_mut(..=N64::new(value)).last();
            if let Some((lower, (lower_prob, ld1, ld2))) = lower {
                let (sign1, param1, sign2, param2) =
                    ProbDistribution1d::calc_sign_params(lp_d2, *lower_prob, prob);
                let (local_d1, local_d2, local_s) = ProbDistribution1d::calc_directions(
                    // lower.raw()-value, lower_prob.ln_sub_exp(prob), upper.raw()-value, upper_prob.ln_sub_exp(prob));
                    // lp_d1-lower.raw(), lp_d2.ln_sub_exp(*lower_prob), value-lower.raw(), prob.ln_sub_exp(*lower_prob));
                    lp_d1 - lower.raw(),
                    param1,
                    sign1,
                    value - lower.raw(),
                    param2,
                    sign2,
                );
                *ld1 = local_s;
                *ld2 = local_d2;
            }
        }

        self.points.insert(N64::new(value), (prob, d1, d2));

        Ok(())
    }

    pub(crate) fn get(&self, value: f64) -> LogProb {
        // println!("size {:?}", self.points.len());
        if self.is_na || value == f64::INFINITY {
            if value == 0.0 {
                // mean or fold change 0
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            let mut upper_it = self.points.range(N64::new(value)..);
            let upper = upper_it.next();
            let upper2 = upper_it.next();
            let mut lower_it = self.points.range(..=N64::new(value)).rev();
            let lower = lower_it.next();
            let lower2 = lower_it.next();
            // println!("value {:?}, lower 2 {:?} upper2  {:?}", value ,lower2, upper2);
            let value = N64::new(value);
            if let Some((upper, (upper_prob, ud1, ud2))) = upper {
                if let Some((lower, (lower_prob, ld1, ld2))) = lower {
                    // Upper and lower bound is there
                    // let factor_u = 1. / (upper.raw() - value).abs();
                    // let factor_l = 1. / (value - lower.raw()).abs();
                    // return (upper_prob + LogProb(factor_u.ln())).ln_add_exp(lower_prob + LogProb(factor_l.ln())) - LogProb((factor_u + factor_l).ln());
                    // println!("value {:?}, lower {:?} upper {:?}", value ,lower, upper);
                    let diff = *upper - *lower;
                    if diff <= 0. {
                        return *upper_prob; // if value is in data structure, upper and lower are the same, no interpolation needed
                    }
                    // println!("diff {:?}", diff);
                    let factor_u = (value - *lower) / diff;
                    let factor_l = (*upper - value) / diff;
                    // println!("factor l {:?} factor u  {:?}", factor_l, factor_u);
                    let min_diff = min((value - lower), (*upper - value));
                    let scaling = (2. * min_diff.raw()) / diff.raw();
                    let scaling_u = ud2 + (LogProb((factor_u.raw() * scaling).ln()));
                    let scaling_l = ld2 + (LogProb((factor_l.raw() * scaling).ln()));
                    let mut result = (LogProb(factor_l.raw().ln()) + lower_prob)
                        .ln_add_exp(LogProb(factor_u.raw().ln()) + upper_prob);

                    if ld1 >= ud1 {
                        //Needed for numerically stable calculations, first add all positive values before doing subtractions
                        if ld1 < &0. {
                            if scaling_l >= result {
                                result = LogProb::ln_zero();
                                return result;
                            }
                            // println!("result {:?} scaling l  {:?}", result, scaling_l);
                            result.ln_sub_exp(scaling_l);
                        } else {
                            result.ln_add_exp(scaling_l);
                        }
                        if ud1 < &0. {
                            if scaling_u >= result {
                                result = LogProb::ln_zero();
                                return result;
                            }
                            // println!("result {:?} scaling u  {:?}", result, scaling_u);
                            result.ln_sub_exp(scaling_u);
                        } else {
                            result.ln_add_exp(scaling_u);
                        }
                    } else {
                        if ud1 < &0. {
                            if scaling_u >= result {
                                result = LogProb::ln_zero();
                                return result;
                            }
                            // println!("result {:?} scaling u  {:?}", result, scaling_u);
                            result.ln_sub_exp(scaling_u);
                        } else {
                            result.ln_add_exp(scaling_u);
                        }
                        if ld1 < &0. {
                            if scaling_l >= result {
                                result = LogProb::ln_zero();
                                return result;
                            }
                            // println!("result {:?} scaling l  {:?}", result, scaling_l);
                            result.ln_sub_exp(scaling_l);
                        } else {
                            result.ln_add_exp(scaling_l);
                        }
                    }
                    //     return (LogProb(factor_l.raw().ln())
                    //             + (lower_prob.ln_add_exp(left))
                    //         .ln_add_exp(LogProb(factor_u.raw().ln())
                    //             + (upper_prob.ln_add_exp(*ud2 + LogProb(scaling.ln())))));
                    return result;
                }
                // Only upper bound
                return upper_prob - LogProb((*upper - value + 1.).raw().ln());
            } else {
                // Only lower bound
                if let Some((lower, (lower_prob, ld1, ld2))) = lower {
                    return lower_prob - LogProb((value - lower + 1.).raw().ln());
                }
                // Neither upper nor lower bound there; should not be reachable
                return LogProb::ln_zero();
            }
        }
    }

    pub(crate) fn normalize(&mut self) -> LogProb {
        if (self.is_na) {
            return LogProb::ln_one();
        }
        // let density = |i, value| self.get(value);
        // let marginals = self
        //     .kdtree
        //     .iter_nearest(&[0.], &squared_euclidean)
        //     .unwrap()
        //     .map(|x| x.0)
        //     .collect::<Vec<_>>()
        //     .windows(2)
        //     .map(|x| LogProb::ln_simpsons_integrate_exp(density, x[0], x[1], 3))
        //     .collect::<Vec<_>>();
        // let marginal = LogProb::ln_sum_exp(&marginals);
        // self.kdtree
        //     .iter_nearest_mut(&[0.], &squared_euclidean)
        //     .unwrap()
        //     .for_each(|x| *x.1 -= marginal );
        // for (prob, d1,d2) in self.points.values_mut() {
        //     *prob = *prob - marginal //Logspace / -> -
        // // }

        let density = |i, value| self.get(value);
        let marginals = self
            .points
            .keys()
            .map(|value| *value)
            .collect::<Vec<_>>()
            .windows(2)
            .map(|x| LogProb::ln_simpsons_integrate_exp(density, x[0].raw(), x[1].raw(), 3))
            .collect::<Vec<_>>();
        let marginal = LogProb::ln_sum_exp(&marginals);
        // self.kdtree
        //     .iter_nearest_mut(&[0.], &squared_euclidean)
        //     .unwrap()
        //     .for_each(|x| *x.1 -= marginal );
        // for (prob, d1,d2) in self.points.values_mut() {
        //     *prob = *prob - marginal //Logspace / -> -
        // }
        let marginal_map = LogProb::ln_trapezoidal_integrate_grid_exp(
            |i, value| self.points.get(&value).unwrap().0,
            &self.points.keys().map(|value| *value).collect::<Vec<_>>(),
        );
        if marginal != LogProb::ln_zero() {
            for (prob, d1, d2) in self.points.values_mut() {
                *prob = *prob - marginal //Logspace / -> -
            }
        }

        // println!("marginals 1 {:?} btreemap  2{:?}", marginal, marginal_map);
        marginal_map
    }
}
