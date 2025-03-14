use anyhow::{Context, Result};
use bio::stats::LogProb;
use getset::Getters;
use itertools_num::linspace;
use serde_derive::{Deserialize, Serialize};
use statrs::distribution::{Continuous, ContinuousCDF, InverseGamma};
use statrs::statistics::Distribution;
use typed_builder::TypedBuilder;

#[derive(TypedBuilder, Copy, Clone, Debug, Getters, Serialize, Deserialize)]
pub(crate) struct PriorParameters {
    shape: f64,
    shift: f64,
    scale: f64,
}

#[derive(Debug, Getters)]
pub(crate) struct Prior {
    inv_gamma: InverseGamma,
    shift: f64,
    #[get = "pub(crate)"]
    left_window: Vec<f64>,
    #[get = "pub(crate)"]
    right_window: Vec<f64>,
}

impl Prior {
    /// Initialize inverse gamma prior. alpha=shape, beta=scale or rate.
    pub(crate) fn new(parameters: &PriorParameters) -> Result<Self> {
        let inv_gamma = InverseGamma::new(parameters.shape, parameters.scale).context(format!(
            "invalid parameters for prior (inverse Gamma) distribution: {:?}",
            parameters
        ))?;
        let window = |a, b, left_exclusive: bool| {
            linspace(a, b, 5)
                .map(|x| x + parameters.shift)
                .skip(if left_exclusive { 1 } else { 0 })
                .collect()
        };
        let left_window = window(
            inv_gamma.inverse_cdf(0.001),
            10., //inv_gamma.mean().unwrap(),
            false,
        );
        let right_window = window(
            10., //inv_gamma.mean().unwrap(),
            20., //inv_gamma.inverse_cdf(0.99), // TODO parameter  mit default 0.99
            true,
        );

        Ok(Prior {
            inv_gamma,
            shift: parameters.shift,
            left_window,
            right_window,
        })
    }

    pub(crate) fn prob(&self, x: f64) -> LogProb {
        LogProb(self.inv_gamma.ln_pdf(x - self.shift)) //TODO
    }

    #[allow(unused)]
    pub(crate) fn mean(&self) -> f64 {
        self.inv_gamma.mean().unwrap() + self.shift
    }

    // pub(crate) fn min_value(&self) -> f64 {
    //     self.inv_gamma.inverse_cdf(0.001)
    //     // self.shift
    // }

    // pub(crate) fn max_value(&self) -> f64 {
    //     20.
    //     //self.inv_gamma.inverse_cdf(0.99)
    // }
}

// #[derive(Debug, Getters)]
// pub(crate) struct Prior {
//     estimated_dispersion: f64,
//     #[get = "pub(crate)"]
//     left_window: Vec<f64>,
//     #[get = "pub(crate)"]
//     right_window: Vec<f64>,
// }

// impl Prior {
//     /// Initialize inverse gamma prior. alpha=shape, beta=scale or rate.
//     pub(crate) fn new(parameter: f64) -> Result<Self> {
//         let window = |a, b, left_exclusive: bool| {
//             linspace(a, b, 5)
//                 // .map(|x| x + parameters.shift)
//                 .skip(if left_exclusive { 1 } else { 0 })
//                 .collect()
//         };
//         let left_window = window(
//             parameter/2.,
//             parameter,
//             false,
//         );
//         let right_window = window(
//             parameter,
//             parameter *2.,
//             true);
//         Ok(Prior {estimated_dispersion:parameter,
//             left_window,
//             right_window,})
//     }

//     pub(crate) fn prob(&self, x: f64) -> LogProb {
//         if x == self.estimated_dispersion {
//             LogProb::ln_one()
//         }else {
//             LogProb::ln_zero()
//         }

//     }

//     #[allow(unused)]
//     pub(crate) fn mean(&self) -> f64 {
//         self.estimated_dispersion
//     }

//     pub(crate) fn min_value(&self) -> f64 {
//         self.estimated_dispersion /2.
//     }

//     pub(crate) fn max_value(&self) -> f64 {
//         self.estimated_dispersion *2.
//     }
// }
