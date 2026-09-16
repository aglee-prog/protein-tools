use crate::model::ProteinResult;
use rusqlite::{Connection, OptionalExtension, params};
use std::{
    path::Path,
    sync::{Arc, Mutex},
    time::{SystemTime, UNIX_EPOCH},
};

#[derive(Clone)]
pub struct Cache {
    connection: Arc<Mutex<Connection>>,
    ttl: i64,
}
impl Cache {
    pub fn open(path: &str, ttl_days: u64) -> Result<Self, String> {
        if let Some(parent) = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
        {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        let connection = Connection::open(path).map_err(|e| e.to_string())?;
        connection
            .busy_timeout(std::time::Duration::from_secs(5))
            .map_err(|e| e.to_string())?;
        connection.execute_batch("PRAGMA journal_mode=WAL; CREATE TABLE IF NOT EXISTS cache (query TEXT NOT NULL, organism TEXT NOT NULL, saved INTEGER NOT NULL, result TEXT NOT NULL, PRIMARY KEY(query, organism));").map_err(|e| e.to_string())?;
        let ttl = ttl_days
            .checked_mul(86400)
            .and_then(|x| i64::try_from(x).ok())
            .ok_or("invalid cache TTL")?;
        Ok(Self {
            connection: Arc::new(Mutex::new(connection)),
            ttl,
        })
    }
    pub async fn get(&self, key: (String, String)) -> Result<Option<ProteinResult>, String> {
        let this = self.clone();
        tokio::task::spawn_blocking(move || {
            let conn = this.connection.lock().map_err(|e| e.to_string())?;
            let value: Option<String> = conn.query_row("SELECT result FROM cache WHERE query=?1 AND organism=?2 AND saved>?3 AND saved<=?4", params![key.0, key.1, now()-this.ttl, now()], |row| row.get(0)).optional().map_err(|e| e.to_string())?;
            value.map(|v| serde_json::from_str(&v).map_err(|e| e.to_string())).transpose()
        }).await.map_err(|e| e.to_string())?
    }
    pub async fn put(&self, key: (String, String), result: ProteinResult) -> Result<(), String> {
        if !result.found || result.error.is_some() {
            return Ok(());
        }
        let this = self.clone();
        tokio::task::spawn_blocking(move || {
            let value = serde_json::to_string(&result).map_err(|e| e.to_string())?;
            this.connection
                .lock()
                .map_err(|e| e.to_string())?
                .execute(
                    "INSERT OR REPLACE INTO cache VALUES (?1, ?2, ?3, ?4)",
                    params![key.0, key.1, now(), value],
                )
                .map_err(|e| e.to_string())?;
            Ok(())
        })
        .await
        .map_err(|e| e.to_string())?
    }
}
fn now() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs() as i64
}
