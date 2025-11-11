use bio::stats::LogProb;
use chrono::offset::Local;
use chrono::DateTime;
use duckdb::{params, AccessMode, Config, Connection};
use itertools::iproduct;
use ordered_float::OrderedFloat;
use std::collections::HashMap;

/// Represents a 2D probability distribution for a given feature,
/// stored in DuckDB.
pub struct ProbDistribution2d {
    conn: Connection,
    feature: String,
    is_na: bool,
}

impl ProbDistribution2d {
    /// Open (and create schema if necessary) for a writable DuckDB connection.
    pub fn open(db_path: &str) -> duckdb::Result<Connection> {
        let conn = Connection::open(db_path)?;
        Self::init_schema(&conn)?;
        Ok(conn)
    }

    /// Ensure the table schema exists.
    pub fn init_schema(conn: &Connection) -> duckdb::Result<()> {
        conn.execute(
            "CREATE TABLE IF NOT EXISTS distributions (
                feature TEXT,
                mu DOUBLE,
                theta DOUBLE,
                prob DOUBLE,
                PRIMARY KEY (feature, mu, theta)
            )",
            [],
        )?;
        Ok(())
    }

    /// Construct a ProbDistribution2d tied to an existing writable connection.
    pub fn with_connection(conn: &Connection, feature: &str) -> duckdb::Result<Self> {
        let is_na = {
            let mut stmt = conn.prepare("SELECT COUNT(*) FROM distributions WHERE feature = ?1")?;
            let mut rows = stmt.query(params![feature])?;
            if let Some(row) = rows.next()? {
                let count: i64 = row.get(0)?;
                count == 0
            } else {
                true
            }
        };
        Ok(Self {
            conn: conn.try_clone()?,
            feature: feature.to_string(),
            is_na,
        })
    }

    /// Open in read-only mode for querying existing results.
    pub fn with_readonly_connection(db_path: &str, feature: &str) -> duckdb::Result<Self> {
        let config = Config::default().access_mode(AccessMode::ReadOnly)?;
        let conn = Connection::open_with_flags(db_path, config)?;
        Ok(Self {
            conn,
            feature: feature.to_string(),
            is_na: false,
        })
    }

    /// Construct an explicit NA distribution (represents missing data).
    pub fn na(feature: &str) -> Self {
        Self {
            conn: Connection::open_in_memory().unwrap(),
            feature: feature.to_string(),
            is_na: true,
        }
    }

    /// Whether this distribution represents a missing feature.
    pub fn is_na(&self) -> bool {
        self.is_na
    }

    /// Compute the full probability grid in memory given mus and thetas.
    pub fn compute_grid<F>(&self, mus: &[f64], thetas: &[f64], mut calc: F) -> Vec<(f64, f64, f64)>
    where
        F: FnMut(f64, f64) -> LogProb,
    {
        let mut results = Vec::with_capacity(mus.len() * thetas.len());
        let total = mus.len() * thetas.len();
        println!(
            "feature {:?} compute_grid total points {}",
            self.feature, total
        );

        let mut count = 0;
        for (j, i) in iproduct!(0..thetas.len(), 0..mus.len()) {
            let mu = mus[i];
            let theta = thetas[j];
            let prob = calc(mu, theta);
            results.push((mu, theta, f64::from(prob)));
            count += 1;

            if count % 10000 == 0 {
                println!(
                    "feature {:?} compute_grid progress {}/{}",
                    self.feature, count, total
                );
            }
        }

        results
    }

    /// Write a computed grid to DuckDB.
    pub fn write_output(&mut self, grid: &[(f64, f64, f64)]) -> duckdb::Result<()> {
        let time1 = std::time::SystemTime::now();
        println!(
            "feature {:?} write_output started at {}",
            self.feature,
            DateTime::<Local>::from(time1).format("%d/%m/%Y %T")
        );

        let tx = self.conn.unchecked_transaction()?;
        let mut stmt = tx.prepare(
            "INSERT INTO distributions (feature, mu, theta, prob)
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT (feature, mu, theta)
             DO UPDATE SET prob = excluded.prob",
        )?;

        for (mu, theta, prob) in grid {
            stmt.execute(params![self.feature, mu, theta, prob])?;
        }

        tx.commit()?;
        let time2 = std::time::SystemTime::now();
        println!(
            "feature {:?} write_output finished at {}, duration {:?}",
            self.feature,
            DateTime::<Local>::from(time2).format("%d/%m/%Y %T"),
            time2.duration_since(time1).unwrap()
        );
        Ok(())
    }

    /// Load a precomputed lookup table for this feature from the DB.
    pub fn load_lookup_table(
        &self,
    ) -> duckdb::Result<HashMap<(OrderedFloat<f64>, OrderedFloat<f64>), LogProb>> {
        let mut stmt = self
            .conn
            .prepare("SELECT mu, theta, prob FROM distributions WHERE feature = ?1")?;
        let mut rows = stmt.query(params![self.feature])?;

        let mut table = HashMap::new();
        while let Some(row) = rows.next()? {
            let mu: f64 = row.get(0)?;
            let theta: f64 = row.get(1)?;
            let prob: f64 = row.get(2)?;
            table.insert((OrderedFloat(mu), OrderedFloat(theta)), LogProb::from(prob));
        }

        Ok(table)
    }

    /// Query a single (mu, theta) pair directly from DuckDB.
    pub fn get(&self, mu: f64, theta: f64) -> LogProb {
        if self.is_na {
            if mu == 0.0 {
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            let mut stmt = self
                .conn
                .prepare(
                    "SELECT prob FROM distributions
                     WHERE feature = ?1 AND mu = ?2 AND theta = ?3",
                )
                .expect("prepare failed");

            let mut rows = stmt.query(params![self.feature, mu, theta]).unwrap();
            if let Some(row) = rows.next().unwrap() {
                let raw_prob: f64 = row.get(0).unwrap();
                LogProb::from(raw_prob)
            } else {
                LogProb::ln_zero()
            }
        }
    }
}
