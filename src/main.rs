use std::path::PathBuf;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use rayon;
use structopt::StructOpt;

mod common;
mod query_points;
mod diff_exp;
mod errors;
mod group_expression;
mod kallisto;
mod preprocess;
mod reduce_features;
mod prior;
mod prob_distribution_1d;
mod prob_distribution_2d;
mod sample_expression;
mod write_fold_changes;

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
            short = "c",
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,
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
        name = "reduce-features",
        about = "Reduce features in datasat to only those contained in list $feature-ids",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    ReduceFeatures {
        #[structopt(
            parse(from_os_str),
            long = "preprocessing_path",
            short = "p",
            help = "Path to preprocessed Kallisto results."
        )]
        preprocessing_path: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "feature-ids",
            short = "i",
            help = "Path to list of feature ids."
        )]
        feature_ids: PathBuf,
        // #[structopt(
        //     parse(from_os_str),
        //     long = "output",
        //     short = "o",
        //     help = "Path to output directory."
        // )]
        // out_dir: PathBuf,
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
            short = "c",
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,
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
            short = "c",
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,
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
            short = "t2",
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
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,
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
        name = "to-text",
        about = "write fold changes from differential expression posteriors between groups into csv file.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    ToText {
        #[structopt(
            parse(from_os_str),
            long = "diff_exp_path",
            short = "d",
            help = "Path to differential expressions"
        )]
        diff_exp_path: PathBuf,
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
            help = "Path to output file."
        )]
        out_file: PathBuf,
    },
    #[structopt(
        name = "kallisto-values",
        about = "write counts or fold changes from kallisto between groups into csv file.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    KallistoValues {
        //add boolean parameter for fold change or counts
        #[structopt(
            long = "foldchange",
            short = "f",
            help = "If --foldchange is set, fold changes are calculated and written into the output file. 
            If --foldchange is not set, counts are written into the output file."
        )]
        foldchange: bool,
        #[structopt(
            parse(from_os_str),
            long = "preprocessing_path",
            short = "p",
            help = "Path to preprocessed Kallisto results."
        )]
        preprocessing_path: PathBuf,
        #[structopt(long = "sample-id", help = "ID of sample to process.")]
        sample_ids: Vec<String>,
        #[structopt(
            parse(from_os_str),
            long = "output",
            short = "o",
            help = "Path to output file."
        )]
        out_file: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::from_args();
    match cli {
        Cli::Preprocess {
            c,
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
            preprocess::preprocess(c, &kallisto_quants, &sample_ids, prior_parameters)
        }
        Cli::ReduceFeatures {
            preprocessing_path,
            feature_ids,
            // out_dir,
        } => {            
            reduce_features::reduce_features(&preprocessing_path, &feature_ids) //, &out_dir)
        }
        Cli::SampleExp {
            preprocessing_path,
            epsilon,
            sample_id,
            c,
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
                c,
                &out_dir,
            )
        }
        Cli::GroupExp {
            preprocessing_path,
            c,
            out_dir,
            threads,
            sample_exprs,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            // calculate per group posteriors
            group_expression::group_expression(&preprocessing_path, &sample_exprs,c, &out_dir)
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
        Cli::ToText {
            diff_exp_path,
            preprocessing_path,
            out_file,
        } => {
            write_fold_changes::write_fold_changes(&preprocessing_path, &diff_exp_path, &out_file)
        }
        Cli::KallistoValues {
            foldchange,
            preprocessing_path,
            sample_ids,
            out_file,
        } => {
            if foldchange {
                write_fold_changes::write_kallisto_fold_changes(&preprocessing_path, sample_ids, &out_file)
            } else {
                write_fold_changes::write_kallisto_counts(&preprocessing_path, sample_ids, &out_file)
            }

        }
    }
}
