use itertools_num::linspace;
use statrs::distribution::{ContinuousCDF, InverseGamma};

pub(crate) struct Prior {
    inv_gamma: InverseGamma,
}

impl Prior {
    pub(crate) fn window(&self) -> impl Iterator<Item = f64> {
        linspace(
            self.inv_gamma.inverse_cdf(0.05),
            self.inv_gamma.inverse_cdf(0.95),
            10,
        ) // TODO how many steps are a good default?
    }
}
