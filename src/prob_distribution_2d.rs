use bio::stats::LogProb;
use duckdb::{params, AccessMode, Config, Connection};
use itertools::iproduct;
use ordered_float::OrderedFloat;
use std::collections::HashMap;
use duckdb::ToSql;

/// Represents a 2D probability distribution for a given feature stored in DuckDB.
pub struct ProbDistribution2d {
    conn: Connection,
    feature: String,
}

pub enum SchemaMode {
    Temp,
    Final,
}

impl ProbDistribution2d {
    // Initialization & Connections

    /// Opens a writable DuckDB connection and ensures schema exists.
    pub fn open(db_path: &str) -> duckdb::Result<Connection> {
        let conn = Connection::open(db_path)?;
        Self::init_schema(&conn)?;
        Ok(conn)
    }

    // /// Ensures that the required schema exists in DuckDB.
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

    /// Attaches a writer to an existing open (writable) connection.
    pub fn with_connection(conn: &Connection, feature: &str) -> duckdb::Result<Self> {
        Ok(Self {
            conn: conn.try_clone()?,
            feature: feature.to_string(),
        })
    }

    /// Opens a DuckDB file in read-only mode for reading stored results.
    pub fn with_readonly_connection(db_path: &str, feature: &str) -> duckdb::Result<Self> {
        let config = Config::default().access_mode(AccessMode::ReadOnly)?;
        let conn = Connection::open_with_flags(db_path, config)?;
        Ok(Self {
            conn,
            feature: feature.to_string(),
        })
    }

    // Writer Thread: Output

    /// Writes a precomputed grid to DuckDB.
    pub fn write_output(&mut self, grid: &[(f64, f64, f64)]) -> duckdb::Result<()> {
        let mut appender = self.conn.appender("distributions")?;

        let feature = &self.feature as &dyn ToSql;

        for (mu, theta, prob) in grid {
            appender.append_row(&[
                feature,
                mu as &dyn ToSql,
                theta as &dyn ToSql,
                prob as &dyn ToSql,
            ])?;
        }
        appender.flush()?;
        Ok(())
    }

    // Reader Thread: Lookup Loading

    /// Loads all stored probabilities for this feature into an in-memory HashMap.
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

    /// Directly query a single (mu, theta) pair from DuckDB.
    pub fn get(&self, mu: f64, theta: f64) -> LogProb {
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

/// Computes a full probability grid in memory, given mus, thetas, and a calc closure.
///
/// This is used by parallel compute threads — no DB connection needed.
pub fn compute_grid<F>(mus: &[f64], thetas: &[f64], mut calc: F) -> Vec<(f64, f64, f64)>
where
    F: FnMut(f64, f64, usize) -> LogProb + Send + Sync,
{
    let mut results = Vec::with_capacity(mus.len() * thetas.len());
    for (j, i) in iproduct!(0..thetas.len(), 0..mus.len()) {
        let mu = mus[i];
        let theta = thetas[j];
        let prob = calc(mu, theta, j);
        results.push((mu, theta, f64::from(prob)));
    }
    results
}
