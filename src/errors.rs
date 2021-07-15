use thiserror::Error;

#[derive(Error, Debug)]
pub(crate) enum Error {
    #[error(
        "Not enough Kallisto quantification files given for normalization (at least 2 required)."
    )]
    NotEnoughQuants,
}
