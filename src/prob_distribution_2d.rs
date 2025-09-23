use duckdb::{Connection, params, Config, AccessMode};
use bio::stats::LogProb;
use itertools::iproduct;
use std::sync::{Arc, Mutex};

/// A probability distribution for a single feature, backed by DuckDB.

pub struct ProbDistribution2d {
    conn: Arc<Mutex<Connection>>, // shared, thread-safe
    feature: String,
    is_na: bool,
}

impl ProbDistribution2d {
    /// Open a shared, mutex-protected connection and initialize schema once.
    pub fn open_shared(db_path: &str) -> duckdb::Result<Arc<Mutex<Connection>>> {
        let conn = Connection::open(db_path)?;
        Self::init_schema(&conn)?; // create table once
        Ok(Arc::new(Mutex::new(conn)))
    }

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

    /// Construct an instance that uses an existing shared connection.
    pub fn with_connection(
        conn: Arc<Mutex<Connection>>,
        feature: &str,
    ) -> duckdb::Result<Self> {
        let is_na = {
            let mut guard = conn.lock().unwrap();
            let mut stmt = guard.prepare("SELECT COUNT(*) FROM distributions WHERE feature = ?1")?;
            let mut rows = stmt.query(params![feature])?;
            if let Some(row) = rows.next()? {
                let count: i64 = row.get(0)?;
                count == 0
            } else {
                true
            }
        };
        Ok(ProbDistribution2d {
            conn,
            feature: feature.to_string(),
            is_na,
        })
    }

    /// Open a DuckDB file in read-only mode for an existing feature
    pub fn with_readonly_connection(db_path: &str, feature: &str) -> duckdb::Result<Self> {
        let config = Config::default().access_mode(AccessMode::ReadOnly)?;
        let conn = Connection::open_with_flags(db_path, config)?;
        Ok(Self {
            conn: Arc::new(Mutex::new(conn)),
            feature: feature.to_string(),
            is_na: false,
        })
    }


   /// Construct an explicit NA distribution (no data for this feature).
    pub fn na(conn: Arc<Mutex<Connection>>, feature: &str) -> Self {
        ProbDistribution2d {
            conn,
            feature: feature.to_string(),
            is_na: true,
        }
    }

    pub fn is_na(&self) -> bool {
        self.is_na
    }

    /// Compute the full probability grid in memory for given mus/thetas.
    pub fn compute_grid<F>(
        &self,
        mus: &[f64],
        thetas: &[f64],
        mut calc: F,
    ) -> Vec<(f64, f64, f64)>
    where
        F: FnMut(f64, f64) -> LogProb,
    {
        let mut results = Vec::with_capacity(mus.len() * thetas.len());
        for (j, i) in iproduct!(0..thetas.len(), 0..mus.len()) {
            let mu = mus[i];
            let theta = thetas[j];
            let prob = calc(mu, theta);
            results.push((mu, theta, f64::from(prob)));
        }
        results
    }

    /// Write a precomputed grid into DuckDB with thread-safe mutex lock.
    pub fn write_output(
        &self,
        grid: &[(f64, f64, f64)],
    ) -> duckdb::Result<()> {
        let mut guard = self.conn.lock().unwrap();
        let tx = guard.unchecked_transaction()?;
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
        Ok(())
    }

    /// Query a probability by exact (mu, theta).
    pub fn get(&self, mu: f64, theta: f64) -> LogProb {
        if self.is_na {
            if mu == 0.0 {
                LogProb::ln_one()
            } else {
                LogProb::ln_zero()
            }
        } else {
            let mut guard = self.conn.lock().unwrap();
            let mut stmt = guard
                .prepare(
                    "SELECT prob FROM distributions WHERE feature = ?1 AND mu = ?2 AND theta = ?3",
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
