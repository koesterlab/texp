use thiserror::Error;
use std::path::PathBuf;

#[derive(Error, Debug)]
pub(crate) enum Error {
    #[error(
        "Not enough Kallisto quantification files given for normalization (at least 2 required)."
    )]
    NotEnoughQuants,
    #[error("Unknown sample id {sample_id}")]
    UnknownSampleId { sample_id: String },
    #[error("Output directory {path} already exists")]
    ExistingOutputDir { path: PathBuf },
}
