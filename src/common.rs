use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::stats::LogProb;
use derefable::Derefable;
use rmp_serde::{Deserializer, Serializer};
use serde::Deserialize as SerdeDeserialize;
use serde::Serialize as SerdeSerialize;
use itertools_num::linspace;
use getset::Getters;
use itertools::iproduct;

use crate::errors::Error;


#[derive(Debug, Getters)]
pub(crate) struct QueryPoints {
    #[get = "pub(crate)"]
    start_points_mu_ik: Vec<f64> ,
    #[get = "pub(crate)"]
    all_mu_ik: Vec<f64> ,
    #[get = "pub(crate)"]
    thetas: Vec<f64>,
    #[get = "pub(crate)"]
    possible_f: Vec<f64>,
}

impl QueryPoints {
    /// Initialize inverse gamma prior. alpha=shape, beta=scale or rate.
    pub(crate) fn new(c: f64) -> Result<Self> {
        let mut start_points_mu_ik: Vec<f64> = linspace(0., 1., 11).collect();
        start_points_mu_ik.extend( linspace(1.5, 100., 198).step_by(10));
        start_points_mu_ik.extend( linspace(101., 500., 400).step_by(10));
        // start_points_mu_ik.extend( linspace(600., 3000., 25));
        // start_points_mu_ik.extend( linspace(4000., 10000., 7));
        println!("len start_points_mu_ik {:?}", start_points_mu_ik.len());

        let mut possible_f: Vec<f64> = linspace(0.05, 5., 100).step_by(2).collect();
        possible_f.extend( linspace(5., 10., 12).step_by(1));
        possible_f.extend( linspace(10.5, 20., 20).step_by(2));
        println!("len possible_f {:?}", possible_f.len());
        
        let mut thetas: Vec<f64> = linspace(0.1, 1., 10).collect();
        // thetas.extend( linspace(1.5, 10., 18).step_by(2));
        // thetas.extend( linspace(11., 165., 155).step_by(10));
        println!("len thetas {:?}", thetas.len());

        let mut all_mu_ik =  start_points_mu_ik.clone();
        for (f, mu_1) in iproduct!(possible_f.clone(), start_points_mu_ik.clone()) {
            let mu_2 = f * (mu_1 + c) - c;
            if mu_2 > 0. {
                all_mu_ik.push(mu_2);
            }            
        }
        all_mu_ik.sort_by(|a, b| a.partial_cmp(b).unwrap()); // TODO NaN -> panic
        all_mu_ik.dedup();
        println!("len all_mu_ik {:?}", all_mu_ik.len());

        Ok(QueryPoints {
            start_points_mu_ik,
            all_mu_ik,
            thetas,
            possible_f,
        })
    }

}

pub(crate) struct Pair {
    pub left: f64,
    pub right: f64,
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct Point {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug)]
pub(crate) struct Square {
    pub top_left: usize,
    pub top_right: usize,
    pub bot_left: usize,
    pub bot_right: usize,
}

pub(crate) fn difference_to_big( // If difference is > than 1/10.000 of maximum probability
    estimated_value: LogProb,
    calculated_value: LogProb,
    prob_dist_max: LogProb
) -> bool {
    // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
    if estimated_value > calculated_value {
        // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
        // println!("ln_sub-max1 {:?}, thresh {:?}, diffToBig {:?}", estimated_value.ln_sub_exp(calculated_value)- prob_dist_max, LogProb::from(0.01_f64.ln()),estimated_value.ln_sub_exp(calculated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()) );
        return estimated_value.ln_sub_exp(calculated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()); 
    }else {
        // println!("est {:?}, calc {:?}, max {:?}", estimated_value, calculated_value, prob_dist_max);
        // println!("ln_sub-max2 {:?}, thresh {:?}, diffToBig {:?}", calculated_value.ln_sub_exp(estimated_value) - prob_dist_max, LogProb::from(0.01_f64.ln()), calculated_value.ln_sub_exp(estimated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()) );
        return calculated_value.ln_sub_exp(estimated_value) - prob_dist_max > LogProb::from(0.01_f64.ln()); 
    }
    // if (estimated_value.exp() - calculated_value.exp()).abs() / prob_dist_max > 0.0001 
}


//--------------------OutDir--------------------

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


