use std::path::PathBuf;

use anyhow::Result;
use structopt::StructOpt;

mod errors;
mod kallisto;
mod normalize;

#[derive(StructOpt, Debug)]
#[structopt(
    name = "bexp",
    about = "Bayesian framework for gene/transcript expression analysis.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
enum Cli {
    #[structopt(
        name = "normalize",
        about = "Calculate scale factors for each given sample by upper quartile normalization.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Normalize {
        #[structopt(
            parse(from_os_str),
            long = "kallisto-quants",
            help = "Paths to Kallisto HDF5 output for each sample."
        )]
        kallisto_quants: Vec<PathBuf>,
        #[structopt(
            long = "sample-ids",
            help = "Sample IDs to use (for each sample given by --kallisto-quants in same order)."
        )]
        sample_ids: Vec<String>,
    },
    #[structopt(
        name = "sample-expression",
        about = "Calculate sample expression likelihoods.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    SampleExp {
        #[structopt(parse(from_os_str), help = "Path to Kallisto HDF5 output for sample.")]
        kallisto_quant: PathBuf,
    },
    #[structopt(
        name = "group-expression",
        about = "Calculate group expression posteriors.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    GroupExp {
        #[structopt(parse(from_os_str))]
        sample_exps: Vec<PathBuf>,
        #[structopt(
            long = "scale-factors",
            short = "s",
            help = "Path to JSON file with per-sample scale factors."
        )]
        scale_factors: PathBuf,
    },
    #[structopt(
        name = "differential-expression",
        about = "Calculate differential expression posteriors between groups.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    DiffExp {
        #[structopt(parse(from_os_str))]
        groups_exps: Vec<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::from_args();
    match cli {
        Cli::Normalize {
            kallisto_quants,
            sample_ids,
        } => {
            // normalize
            normalize::normalize(&kallisto_quants, &sample_ids)
        }
        Cli::SampleExp { kallisto_quant } => {
            // calculate per sample likelihoods
            Ok(())
        }
        Cli::GroupExp {
            sample_exps,
            scale_factors,
        } => {
            // calculate per group posteriors
            Ok(())
        }
        Cli::DiffExp { groups_exps } => {
            // calculate between group differential expressions
            Ok(())
        }
    }
}
