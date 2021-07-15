use std::path::Path;

use anyhow::Result;
use hdf5;
use ndarray::{s, Array1, Array2, Dim};
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

    pub(crate) fn bootstrapped_counts(&self) -> Result<Array2<f64>> {
        let num_bootstraps = self
            .reader
            .dataset("aux/num_bootstrap")?
            .read_1d::<usize>()?[0];
        let seq_length = self.reader.dataset("aux/lengths")?.read_1d::<f64>()?;

        let mut bootstraps = Array2::zeros(Dim([num_bootstraps, seq_length.len()]));

        for i in 0..num_bootstraps {
            let dataset = self.reader.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let est_counts = dataset.read_1d::<f64>()?;
            let norm_counts = est_counts / &seq_length;
            bootstraps.slice_mut(s![i, ..]).assign(&norm_counts);
        }

        Ok(bootstraps)
    }

    pub(crate) fn feature_ids(&self) -> Result<Array1<String>> {
        Ok(self
            .reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?
            .mapv(|id| id.to_string()))
    }
}
