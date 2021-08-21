use anyhow::Result;
use bio::stats::{LogProb, Prob};
use itertools_num::linspace;
use statrs::distribution::{Continuous, ContinuousCDF, InverseGamma};

pub(crate) struct Prior {
    inv_gamma: InverseGamma,
    shift: f64,
}

impl Prior {
    /// Initialize inverse gamma prior. alpha=shape, beta=scale or rate.
    pub(crate) fn new(shape: f64, scale: f64, shift: f64) -> Result<Self> {
        Ok(Prior {
            inv_gamma: InverseGamma::new(shape, scale)?,
            shift,
        })
    }

    pub(crate) fn prob(&self, x: f64) -> LogProb {
        LogProb::from(Prob(self.inv_gamma.pdf(x - self.shift)))
    }

    pub(crate) fn window(&self) -> impl Iterator<Item = f64> {
        linspace(
            self.inv_gamma.inverse_cdf(0.05),
            self.inv_gamma.inverse_cdf(0.95),
            10,
        ) // TODO how many steps are a good default?
    }
}
