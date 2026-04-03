use bio::stats::LogProb;
use duckdb::{AccessMode, Config, Connection, params};
use ordered_float::OrderedFloat;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
extern crate chrono;

pub struct ProbDistribution1d {
    conn: Arc<Mutex<Connection>>,
    feature: String,
    is_na: bool,
}

impl ProbDistribution1d {
    /// Open a shared, mutex-protected connection and initialize schema once.
    pub fn open_shared(db_path: &str) -> duckdb::Result<Arc<Mutex<Connection>>> {
        let conn = Connection::open(db_path)?;
        Self::init_schema(&conn)?; // create table once
        Ok(Arc::new(Mutex::new(conn)))
    }

    /// Initialize the table schema (call once per DB file)
    pub fn init_schema(conn: &Connection) -> duckdb::Result<()> {
        conn.execute(
            "CREATE TABLE IF NOT EXISTS distributions_1d (
                feature TEXT,
                x DOUBLE,
                prob DOUBLE,
                PRIMARY KEY (feature, x)
            )",
            [],
        )?;
        Ok(())
    }

    /// Construct an instance that uses an existing shared connection.
    pub fn with_connection(conn: Arc<Mutex<Connection>>, feature: &str) -> duckdb::Result<Self> {
        let is_na = {
            let guard = conn.lock().unwrap();
            let mut stmt =
                guard.prepare("SELECT COUNT(*) FROM distributions_1d WHERE feature = ?1")?;
            let mut rows = stmt.query(params![feature])?;
            if let Some(row) = rows.next()? {
                let count: i64 = row.get(0)?;
                count == 0
            } else {
                true
            }
        };

        Ok(Self {
            conn,
            feature: feature.to_string(),
            is_na,
        })
    }

    /// Open a DuckDB file in read-only mode for an existing feature.
    pub fn with_readonly_connection(db_path: &str, feature: &str) -> duckdb::Result<Self> {
        let config = Config::default().access_mode(AccessMode::ReadOnly)?;
        let conn = Connection::open_with_flags(db_path, config)?;
        Ok(Self {
            conn: Arc::new(Mutex::new(conn)),
            feature: feature.to_string(),
            is_na: false,
        })
    }

    /// Insert a single point
    pub fn insert(&self, x: f64, prob: LogProb) -> duckdb::Result<()> {
        let guard = self.conn.lock().unwrap();
        let tx = guard.unchecked_transaction()?;
        let mut stmt = tx.prepare(
            "INSERT INTO distributions_1d (feature, x, prob)
             VALUES (?1, ?2, ?3)
             ON CONFLICT (feature, x) DO UPDATE SET prob = excluded.prob",
        )?;
        stmt.execute(params![self.feature, x, f64::from(prob)])?;
        tx.commit()?;
        Ok(())
    }

    // /// Query a single probability
    // pub fn get(&self, x: f64) -> LogProb {
    //     let mut guard = self.conn.lock().unwrap();
    //     let mut stmt = guard
    //         .prepare("SELECT prob FROM distributions_1d WHERE feature = ?1 AND x = ?2")
    //         .unwrap();
    //     let mut rows = stmt.query(params![self.feature, x]).unwrap();
    //     if let Some(row) = rows.next().unwrap() {
    //         LogProb::from(row.get::<_, f64>(0).unwrap())
    //     } else {
    //         LogProb::ln_zero()
    //     }
    // }
    // /// Insert a single point.
    // pub fn insert(&self, x: f64, prob: LogProb) -> duckdb::Result<()> {
    //     self.insert_points(&[(x, prob)])
    // }

    /// Query a probability for a given x (direct DB query).
    pub fn get(&self, x: f64) -> LogProb {
        let guard = self.conn.lock().unwrap();
        let mut stmt = guard
            .prepare("SELECT prob FROM distributions_1d WHERE feature = ?1 AND x = ?2")
            .unwrap();
        let mut rows = stmt.query(params![self.feature, x]).unwrap();
        if let Some(row) = rows.next().unwrap() {
            let raw_prob: f64 = row.get(0).unwrap();
            LogProb::from(raw_prob)
        } else {
            LogProb::ln_zero()
        }
    }

    /// Load the whole feature into a lookup table for fast repeated queries.
    pub fn load_lookup_table(&self) -> duckdb::Result<HashMap<OrderedFloat<f64>, LogProb>> {
        let guard = self.conn.lock().unwrap();
        let mut stmt = guard.prepare("SELECT x, prob FROM distributions_1d WHERE feature = ?1")?;
        let mut rows = stmt.query(params![self.feature])?;
        let mut out = HashMap::new();

        while let Some(row) = rows.next()? {
            let x: f64 = row.get(0)?;
            let prob: f64 = row.get(1)?;
            // use rounded integer key for hashing (no f64 Eq/Hash issues)
            // let x_key = (x.to_bits() as i64);
            out.insert(OrderedFloat(x), LogProb::from(prob));
        }
        Ok(out)
    }

    //  /// Write a precomputed 1D grid into DuckDB with thread-safe mutex lock
    // pub fn write_output(&self, grid: &[(f64, LogProb)]) -> Result<()> {
    //     let time1 = std::time::SystemTime::now()
    //     println!(
    //         "feature {:?} write_output started at {}",
    //         self.feature,
    //         chrono::DateTime::<Local>::from(time1).format("%d/%m/%Y %T")
    //     );

    //     let mut guard = self.conn.lock().unwrap();
    //     let tx = guard.unchecked_transaction()?;
    //     let mut stmt = tx.prepare(
    //         "INSERT INTO distributions_1d (feature, x, prob)
    //          VALUES (?1, ?2, ?3)
    //          ON CONFLICT (feature, x)
    //          DO UPDATE SET prob = excluded.prob",
    //     )?;

    //     for (x, prob) in grid {
    //         stmt.execute(params![self.feature, *x, f64::from(*prob)])?;
    //     }

    //     tx.commit()?;

    //     let time2 = std::time::SystemTime::now()
    //     println!(
    //         "feature {:?} write_output finished at {}, duration {:?}",
    //         self.feature,
    //         chrono::DateTime::<Local>::from(time2).format("%d/%m/%Y %T"),
    //         time2.duration_since(time1).unwrap()
    //     );

    //     Ok(())
    // }
}
