use anyhow::Result;
use bio::stats::{LogProb, Prob};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

mod common;
mod diff_exp;
mod errors;
mod group_expression;
mod kallisto;
mod preprocess;
mod prior;
mod prob_distribution_1d;
mod prob_distribution_2d;
mod query_points;
mod reduce_features;
mod sample_expression;
mod write_fold_changes;

#[derive(Parser, Debug)]
#[command(
    name = "t-exp",
    about = "Tyrannosaurus Exp: Bayesian framework for gene/transcript expression analysis.",
    version
)]
enum Cli {
    /// Calculate mean and dispersion estimates as well as scale factors for each given sample by upper quartile normalization.
    Preprocess {
        #[arg(
            short,
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,

        #[arg(long, help = "Paths to Kallisto HDF5 output for each sample.")]
        kallisto_quants: Vec<PathBuf>,

        #[arg(
            long,
            help = "Sample IDs to use (for each sample given by --kallisto-quants in same order)."
        )]
        sample_ids: Vec<String>,

        #[arg(
            long,
            default_value = "1.0409428761583088",
            help = "Shape of prior distribution (inverse gamma)."
        )]
        prior_shape: f64,

        #[arg(
            long,
            default_value = "2.064553353135377",
            help = "Scale of prior distribution (inverse gamma)."
        )]
        prior_scale: f64,

        #[arg(
            long,
            default_value = "-0.017934198042149123",
            help = "Shift of prior distribution (inverse gamma)."
        )]
        prior_shift: f64,
    },

    /// Reduce features in dataset to only those contained in list $feature-ids
    ReduceFeatures {
        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(long, short = 'i', help = "Path to list of feature ids.")]
        feature_ids: PathBuf,
    },

    /// Calculate sample expression likelihoods.
    SampleExp {
        #[arg(long, help = "ID of sample to process.")]
        sample_id: String,

        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(long, short = 'o', help = "Path to output directory.")]
        out_dir: PathBuf,

        #[arg(
            long,
            default_value = "1e-9",
            help = "Epsilon for stopping likelihood calculation."
        )]
        epsilon: f64,

        #[arg(
            short,
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,

        #[arg(long, default_value = "1", help = "Number of threads to use.")]
        threads: usize,
    },

    /// Calculate group expression posteriors.
    GroupExp {
        #[arg(help = "Paths to sample expressions.")]
        sample_exprs: Vec<PathBuf>,

        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(
            short,
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,

        #[arg(long, short = 'o', help = "Path to output directory.")]
        output: PathBuf,

        #[arg(long, default_value = "1", help = "Number of threads to use.")]
        threads: usize,
    },

    /// Calculate differential expression posteriors between groups.
    #[command(name = "differential-expression")]
    DiffExp {
        #[arg(
            long = "group_path1",
            short = '1',
            help = "Path to group expressions of group 1."
        )]
        group_path1: PathBuf,

        #[arg(
            long = "group_path2",
            short = '2',
            help = "Path to group expressions of group 2."
        )]
        group_path2: PathBuf,

        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(
            short,
            default_value = "0",
            help = "Pseudo counts c for fold change calculation."
        )]
        c: f64,

        #[arg(long, short = 'o', help = "Path to output directory.")]
        output: PathBuf,

        #[arg(long, default_value = "1", help = "Number of threads to use.")]
        threads: usize,
    },

    /// Write fold changes from differential expression posteriors between groups into csv file.
    #[command(name = "to-text")]
    ToText {
        #[arg(long, short = 'd', help = "Path to differential expressions")]
        diff_exp_path: PathBuf,

        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(long = "output_dist", short = 'd', help = "Path to output file.")]
        output_dist: PathBuf,

        #[arg(
            long = "output_max_prob_fc",
            short = 'm',
            help = "Path to output file."
        )]
        out_file_max_prob_fc: PathBuf,
    },

    /// Write counts or fold changes from kallisto between groups into csv file.
    #[command(name = "kallisto-values")]
    KallistoValues {
        #[arg(
            long,
            short = 'f',
            help = "If --foldchange is set, fold changes are calculated. Otherwise, counts are written."
        )]
        foldchange: bool,

        #[arg(long, short = 'p', help = "Path to preprocessed Kallisto results.")]
        preprocessing_path: PathBuf,

        #[arg(long = "sample-id", help = "ID of sample to process.")]
        sample_ids: Vec<String>,

        #[arg(long, short = 'o', help = "Path to output file.")]
        out_file: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

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
        } => reduce_features::reduce_features(&preprocessing_path, &feature_ids),
        Cli::SampleExp {
            preprocessing_path,
            epsilon,
            sample_id,
            c,
            out_dir,
            threads,
        } => sample_expression::sample_expression(
            &preprocessing_path,
            &sample_id,
            LogProb::from(Prob::checked(epsilon)?),
            c,
            threads,
            &out_dir,
        ),
        Cli::GroupExp {
            preprocessing_path,
            c,
            output,
            threads,
            sample_exprs,
        } => group_expression::group_expression(
            &preprocessing_path,
            &sample_exprs,
            c,
            threads,
            &output,
        ),
        Cli::DiffExp {
            group_path1,
            group_path2,
            preprocessing_path,
            c,
            output: out_dir,
            threads: _, // threads unused in original match arm
        } => diff_exp::diff_exp(c, &preprocessing_path, &group_path1, &group_path2, &out_dir),
        Cli::ToText {
            diff_exp_path,
            preprocessing_path,
            output_dist,
            out_file_max_prob_fc,
        } => write_fold_changes::write_fold_changes(
            &preprocessing_path,
            &diff_exp_path,
            &output_dist,
            &out_file_max_prob_fc,
        ),
        Cli::KallistoValues {
            foldchange,
            preprocessing_path,
            sample_ids,
            out_file,
        } => {
            if foldchange {
                write_fold_changes::write_kallisto_fold_changes(
                    &preprocessing_path,
                    sample_ids,
                    &out_file,
                )
            } else {
                write_fold_changes::write_kallisto_counts(
                    &preprocessing_path,
                    sample_ids,
                    &out_file,
                )
            }
        }
    }
}
