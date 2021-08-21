use anyhow::Result;
use bio::stats::{LogProb, Prob};
use getset::Getters;
use itertools_num::linspace;
use statrs::distribution::{Continuous, ContinuousCDF, InverseGamma};
use statrs::statistics::Distribution;

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
    pub(crate) fn new(shape: f64, scale: f64, shift: f64) -> Result<Self> {
        let inv_gamma = InverseGamma::new(shape, scale)?;
        let window = |a, b, left_exclusive: bool| {
            linspace(a, b, 10)
                .map(|x| x + shift)
                .skip(if left_exclusive { 1 } else { 0 })
                .collect()
        };
        let left_window = window(
            inv_gamma.inverse_cdf(0.05),
            inv_gamma.mean().unwrap(),
            false,
        );
        let right_window = window(inv_gamma.mean().unwrap(), inv_gamma.inverse_cdf(0.95), true);

        Ok(Prior {
            inv_gamma,
            shift,
            left_window,
            right_window,
        })
    }

    pub(crate) fn prob(&self, x: f64) -> LogProb {
        LogProb(self.inv_gamma.ln_pdf(x - self.shift))
    }

    pub(crate) fn mean(&self) -> f64 {
        self.inv_gamma.mean().unwrap() + self.shift
    }

    pub(crate) fn min_value(&self) -> f64 {
        self.shift
    }

    pub(crate) fn max_value(&self) -> f64 {
        self.inv_gamma.inverse_cdf(0.9999)
    }
}
