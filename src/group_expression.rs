//! This implements formula 5, 6, 7 of the document.
use anyhow::Result;
use bio::stats::LogProb;
use duckdb::Connection;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;
use std::sync::mpsc::{Receiver, Sender};
use std::thread;

use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;

pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    c: f64,
    out_dir_path: &Path,
) -> Result<()> {
    // Load preprocessing
    let preprocessing = Preprocessing::from_path(preprocessing)?;
    let sample_ids = preprocessing
        .scale_factors()
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    let feature_ids: Vec<_> = preprocessing.feature_ids().iter().enumerate().collect();

    // Shared writer connection (only writer thread touches this)
    let db_path = out_dir_path.to_str().unwrap().to_string();

    // --- Channel setup ---
    let (tx, rx): (
        Sender<(String, Vec<(f64, f64, f64)>)>,
        Receiver<(String, Vec<(f64, f64, f64)>)>,
    ) = channel();

    // --- Spawn the writer thread ---
    let writer_handle = thread::spawn(move || {
        let conn = Connection::open(&db_path).expect("Failed to open DuckDB in writer");
        ProbDistribution2d::init_schema(&conn).unwrap();

        // Writer loop
        while let Ok((feature_id, probs)) = rx.recv() {
            let mut likelihoods = ProbDistribution2d::with_connection(&conn, &feature_id).unwrap();
            likelihoods.write_output(&probs).unwrap();
        }
    });

    // -----------------------------------------
    // Parallel worker threads (compute only)
    // -----------------------------------------
    feature_ids
        .par_iter()
        .try_for_each(|(i, feature_id)| -> Result<()> {
            println!("Feature {i} ({feature_id}) — starting computation");

            // Open per-sample likelihood tables (read-only)
            let sample_expression_likelihoods: Vec<_> = sample_expression_paths
                .iter()
                .map(|path| {
                    ProbDistribution2d::with_readonly_connection(path.to_str().unwrap(), feature_id)
                })
                .collect::<duckdb::Result<_>>()?;

            // Preload all lookup tables in memory
            let lookup_tables: Vec<_> = sample_expression_likelihoods
                .iter()
                .map(|likelihood| likelihood.load_lookup_table())
                .collect::<duckdb::Result<_>>()?;

            // let prob_dist = ProbDistribution2d::new(feature_id);
            let prob_dist = ProbDistribution2d::na(feature_id);

            let calc_prob = |mu_ik: f64, theta_i: f64| {
                if mu_ik == 0.0 {
                    return LogProb::ln_zero();
                }

                let key = (OrderedFloat(mu_ik), OrderedFloat(theta_i));
                let probs: Vec<LogProb> = lookup_tables
                    .iter()
                    .map(|table| table.get(&key).cloned().unwrap_or(LogProb::ln_zero()))
                    .collect();

                LogProb::ln_sum_exp(&probs)
            };

            let query_points = query_points::calc_query_points(
                c,
                preprocessing.mean_disp_estimates().clone(),
                sample_ids.clone(),
                preprocessing.feature_ids().clone(),
                *i,
            );

            let mu_ik_points = query_points.all_mu_ik();
            let theta_points = query_points.thetas();

            let probs = prob_dist.compute_grid(&mu_ik_points, &theta_points, calc_prob);

            // Send computed result to writer thread
            tx.send((feature_id.to_string(), probs)).unwrap();

            println!("Feature {i} ({feature_id}) — computation done, sent to writer.");
            Ok(())
        })?;

    drop(tx); // close channel
    writer_handle.join().unwrap();

    Ok(())
}
