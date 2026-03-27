//! This implements formula 5, 6, 7 of the document.
use anyhow::Result;
use bio::stats::LogProb;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use rayon::ThreadPoolBuilder;
use std::fs;
use std::collections::HashMap;

use crate::preprocess::Preprocessing;
use crate::prob_distribution_2d::compute_grid;
use crate::prob_distribution_2d::ProbDistribution2d;
use crate::query_points;

pub(crate) fn group_expression(
    preprocessing: &Path,
    sample_expression_paths: &[PathBuf],
    c: f64,
    threads: usize,
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



    let base_db_path = out_dir_path.to_str().unwrap().to_string();

    let query_points_per_feature = query_points::QueryPointsPerFeature::new(&preprocessing, c);


    // -------------------------------
    // 1. HPC SAFEGUARD
    // -------------------------------
    let max_threads = std::cmp::min(threads, 32);
    let custom_pool = ThreadPoolBuilder::new()
        .num_threads(max_threads)
        .build()?;

    // -------------------------------
    // 2. SCATTER PHASE
    // -------------------------------
    let temp_paths: Vec<String> = (0..max_threads)
        .map(|i| format!("{}_temp_{}.duckdb", base_db_path, i))
        .collect();

    custom_pool.install(|| {
        feature_ids
            .into_par_iter()
            // .try_for_each_init(
            .map_init(
                 // INIT → runs once per thread
            {
                let temp_paths = temp_paths.clone();

                move || {
                    let thread_idx = rayon::current_thread_index().unwrap();

                    let path = &temp_paths[thread_idx];

                    let conn = ProbDistribution2d::open(path)
                        .expect("Failed to open DuckDB");

                    (path.clone(), conn)
                }
            },
             // WORK
            |(temp_path, conn), (i, feature_id)| {
                // Open per-sample likelihood tables (read-only)
                let sample_expression_likelihoods: Vec<_> = sample_expression_paths
                    .iter()
                    .map(|path| {
                        ProbDistribution2d::with_readonly_connection(
                            path.to_str().unwrap(),
                            &feature_id,
                        )
                    })
                    .collect::<duckdb::Result<_>>()
                    .unwrap();

                // Preload all lookup tables in memory
                let lookup_tables: Vec<_> = sample_expression_likelihoods
                    .iter()
                    .map(|likelihood| match likelihood.load_lookup_table() {
                        Ok(table) => table,
                        Err(e) => {
                            eprintln!("Failed to load lookup table: {:?}", e);
                            HashMap::new() // or whatever fallback
                        }
                    })
                    .collect();

                let calc_prob = |mu_ik: f64, theta_i: f64, theta_idx: usize| {
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

                let query_points = query_points_per_feature.get(i);

                let mu_ik_points = query_points.all_mu_ik();
                let theta_points = query_points.thetas();

                let probs = compute_grid(&mu_ik_points, &theta_points, calc_prob);

                let mut writer = ProbDistribution2d::with_connection(&*conn, &feature_id.to_string()).expect("writer");

                writer.write_output(&probs).expect("write failed");
            },
        )
        .count(); // Force execution of the parallel iterator
    });


    // -------------------------------
    // 3. GATHER PHASE
    // -------------------------------
    let final_conn = duckdb::Connection::open(&base_db_path)?;
    ProbDistribution2d::init_schema(&final_conn)?;

    for temp_path in &temp_paths {
        final_conn.execute(
            &format!("ATTACH '{}' AS temp_db", temp_path),
            [],
        )?;

        final_conn.execute(
            "INSERT OR REPLACE INTO distributions SELECT * FROM temp_db.distributions",
            [],
        )?;

        final_conn.execute("DETACH temp_db", [])?;
    }
    // cleanup phase (parallel)
    temp_paths.par_iter().for_each(|path| {
        let _ = fs::remove_file(path);
    });

    Ok(())
}
