use std::path::Path;

use anyhow::Result;
use hdf5;
use ndarray::Array1;
use noisy_float::types::N64;

pub(crate) struct KallistoQuant {
    reader: hdf5::File,
}

impl KallistoQuant {
    pub(crate) fn new(kallisto_quant: &Path) -> Result<Self> {
        Ok(KallistoQuant {
            reader: hdf5::File::open(kallisto_quant)?,
        })
    }

    pub(crate) fn len_norm_counts(&self) -> Result<Array1<N64>> {
        let counts = self.reader.dataset("est_counts")?.read_1d::<f64>()?;
        let lens = self.reader.dataset("aux/lengths")?.read_1d::<f64>()?;
        Ok((counts / lens).mapv(|v| N64::unchecked_new(v)))
    }
}
