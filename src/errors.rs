use std::path::PathBuf;
use thiserror::Error;

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
    #[error("Input directory {path} does not exist")]
    NotExistingOutputDir { path: PathBuf },
}
