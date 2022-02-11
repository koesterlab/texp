use std::path::PathBuf;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use rayon;
use structopt::StructOpt;

use crate::common::{MeanDispersionPair, Outdir};

mod common;
mod diff_exp;
mod errors;
mod group_expression;
mod kallisto;
mod preprocess;
mod prior;
mod sample_expression;
use common::ProbDistribution;

#[derive(StructOpt, Debug)]
#[structopt(
    name = "t-exp",
    about = "Tyrannosaurus Exp: Bayesian framework for gene/transcript expression analysis.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
enum Cli {
    #[structopt(
        name = "preprocess",
        about = "Calculate mean and dispersion estimates as well as scale factors for each given sample by upper quartile normalization.",
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
        #[structopt(
            parse(from_os_str),
            long = "preprocessing_path",
            short = "p",
            help = "Path to preprocessed Kallisto results."
        )]
        preprocessing_path: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "output",
            short = "o",
            help = "Path to output directory."
        )]
        out_dir: PathBuf,
        #[structopt(
            long = "epsilon",
            default_value = "1e-9",
            help = "Epsilon for stopping likelihood calculation."
        )]
        epsilon: f64,
        #[structopt(
            long = "threads",
            default_value = "1",
            help = "Number of threads to use."
        )]
        threads: usize,
    },
    #[structopt(
        name = "group-expression",
        about = "Calculate group expression posteriors.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    GroupExp {
        #[structopt(parse(from_os_str), help = "Paths to sample expressions.")]
        sample_exprs: Vec<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long = "preprocessing_path",
            short = "p",
            help = "Path to preprocessed Kallisto results."
        )]
        preprocessing_path: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "output",
            short = "o",
            help = "Path to output directory."
        )]
        out_dir: PathBuf,
        #[structopt(
            long = "threads",
            default_value = "1",
            help = "Number of threads to use."
        )]
        threads: usize,
    },
    #[structopt(
        name = "differential-expression",
        about = "Calculate differential expression posteriors between groups.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    DiffExp {
        #[structopt(
            parse(from_os_str),
            long = "group_path1",
            short = "g1",
            help = "Path to group expressions of group 1."
        )]
        group_path1: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "group_path2",
            short = "g2",
            help = "Path to group expressions of group 2."
        )]
        group_path2: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "preprocessing_path",
            short = "p",
            help = "Path to preprocessed Kallisto results."
        )]
        preprocessing_path: PathBuf,
        #[structopt(
            short = "c", 
            default_value = "10", 
            help = "Pseudo counts c for fold change calculation.")]
        c: f32,
        #[structopt(
            parse(from_os_str),
            long = "output",
            short = "o",
            help = "Path to output directory."
        )]
        out_dir: PathBuf,
        #[structopt(
            long = "threads",
            default_value = "1",
            help = "Number of threads to use."
        )]
        threads: usize,
    },
    #[structopt(
        name = "show-sample-expression",
        about = "Decode mpk into JSON.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    ShowSampleExpressions {
        #[structopt(parse(from_os_str), help = "Path to sample expressions dir.")]
        path: PathBuf,
        #[structopt(long = "feature-id", help = "ID of feature to show.")]
        feature_id: String,
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
            let prior_parameters = prior::PriorParameters::builder()
                .shape(prior_shape)
                .scale(prior_scale)
                .shift(prior_shift)
                .build();
            // normalize
            preprocess::preprocess(&kallisto_quants, &sample_ids, prior_parameters)
        }
        Cli::SampleExp {
            preprocessing_path,
            epsilon,
            sample_id,
            out_dir,
            threads,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            // calculate per sample likelihoods
            sample_expression::sample_expression(
                &preprocessing_path,
                &sample_id,
                LogProb::from(Prob::checked(epsilon)?),
                &out_dir,
            )
        }
        Cli::GroupExp {
            preprocessing_path,
            out_dir,
            threads,
            sample_exprs,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            // calculate per group posteriors
            group_expression::group_expression(&preprocessing_path, &sample_exprs, &out_dir)
        }
        Cli::DiffExp {
            group_path1,
            group_path2,
            preprocessing_path,
            c,
            out_dir,
            threads,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
            // calculate differential expression between groups
            diff_exp::diff_exp(c, &preprocessing_path, &group_path1, &group_path2, &out_dir)
        }
        Cli::ShowSampleExpressions { path, feature_id } => {
            let dir = Outdir::open(&path)?;
            let expr: ProbDistribution<MeanDispersionPair> = dir.deserialize_value(&feature_id)?;
            // TODO output as tsv (use csv crate)
            dbg!(&expr);
            Ok(())
        }
    }
}
