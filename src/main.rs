use std::path::PathBuf;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use structopt::StructOpt;

mod common;
mod errors;
mod group_expression;
mod kallisto;
mod preprocess;
mod prior;
mod sample_expression;

#[derive(StructOpt, Debug)]
#[structopt(
    name = "bexp",
    about = "Bayesian framework for gene/transcript expression analysis.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
enum Cli {
    #[structopt(
        name = "preprocess",
        about = "Calculate scale factors for each given sample by upper quartile normalization.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Preprocess {
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
        #[structopt(
            long = "prior-shape",
            default_value = "1.0409428761583088",
            help = "Shape of prior distribution (inverse gamma)."
        )]
        prior_shape: f64,
        #[structopt(
            long = "prior-scale",
            default_value = "2.064553353135377",
            help = "Scale of prior distribution (inverse gamma)."
        )]
        prior_scale: f64,
        #[structopt(
            long = "prior-shift",
            default_value = "-0.017934198042149123",
            help = "Shift of prior distribution (inverse gamma)."
        )]
        prior_shift: f64,
    },
    #[structopt(
        name = "sample-expression",
        about = "Calculate sample expression likelihoods.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    SampleExp {
        #[structopt(long = "sample-id", help = "ID of sample to process.")]
        sample_id: String,
        #[structopt(parse(from_os_str), help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,
        #[structopt(
            long = "epsilon",
            default_value = "1e-9",
            help = "Epsilon for stopping likelihood calculation."
        )]
        epsilon: f64,
    },
    #[structopt(
        name = "group-expression",
        about = "Calculate group expression posteriors.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    GroupExp {
        #[structopt(parse(from_os_str), help = "Paths to sample expressions.")]
        sample_exprs: Vec<PathBuf>,
        #[structopt(parse(from_os_str), help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,
        #[structopt(parse(from_os_str), help = "Path to output directory.")]
        out_dir: PathBuf,
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
        Cli::Preprocess {
            kallisto_quants,
            sample_ids,
            prior_shape,
            prior_scale,
            prior_shift,
        } => {
            let prior_parameters =
                prior::PriorParameters::new(prior_shape, prior_scale, prior_shift);
            // normalize
            preprocess::preprocess(&kallisto_quants, &sample_ids, prior_parameters)
        }
        Cli::SampleExp {
            preprocessing_path,
            epsilon,
            sample_id,
        } => {
            // calculate per sample likelihoods
            sample_expression::sample_expression(
                &preprocessing_path,
                &sample_id,
                LogProb::from(Prob::checked(epsilon)?),
            )
        }
        Cli::GroupExp {
            sample_exprs,
            preprocessing_path,
            out_dir,
        } => {
            // calculate per group posteriors
            group_expression::group_expression(&preprocessing_path, &sample_exprs, &out_dir)
        }
        Cli::DiffExp { groups_exps } => {
            // calculate between group differential expressions
            Ok(())
        }
    }
}
