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
        connection.execute_batch("PRAGMA journal_mode=WAL; CREATE TABLE IF NOT EXISTS cache (query TEXT NOT NULL, organism TEXT NOT NULL, saved INTEGER NOT NULL, result TEXT NOT NULL, PRIMARY KEY(query, organism)); CREATE TABLE IF NOT EXISTS kegg_cache (organism TEXT NOT NULL, stage TEXT NOT NULL, identifier TEXT NOT NULL, saved INTEGER NOT NULL, result TEXT NOT NULL, PRIMARY KEY(organism, stage, identifier));").map_err(|e| e.to_string())?;
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
    pub async fn get_kegg(
        &self,
        organism: &str,
        stage: &str,
        identifier: &str,
    ) -> Result<Option<Vec<String>>, String> {
        self.get_kegg_value(organism, stage, identifier).await
    }
    pub async fn put_kegg(
        &self,
        organism: &str,
        stage: &str,
        identifier: &str,
        values: Vec<String>,
    ) -> Result<(), String> {
        self.put_kegg_value(organism, stage, identifier, values)
            .await
    }
    // KEGG stages store normalized values per identifier, including successful empty mappings.
    pub async fn get_kegg_value<T: serde::de::DeserializeOwned + Send + 'static>(
        &self,
        organism: &str,
        stage: &str,
        identifier: &str,
    ) -> Result<Option<T>, String> {
        let this = self.clone();
        let key = (organism.to_owned(), stage.to_owned(), identifier.to_owned());
        tokio::task::spawn_blocking(move || {
            let conn = this.connection.lock().map_err(|e| e.to_string())?;
            let value: Option<String> = conn.query_row(
                "SELECT result FROM kegg_cache WHERE organism=?1 AND stage=?2 AND identifier=?3 AND saved>?4 AND saved<=?5",
                params![key.0, key.1, key.2, now()-this.ttl, now()], |row| row.get(0)
            ).optional().map_err(|e| e.to_string())?;
            value.map(|v| serde_json::from_str(&v).map_err(|e| e.to_string())).transpose()
        }).await.map_err(|e| e.to_string())?
    }
    pub async fn put_kegg_value<T: serde::Serialize + Send + 'static>(
        &self,
        organism: &str,
        stage: &str,
        identifier: &str,
        values: T,
    ) -> Result<(), String> {
        let this = self.clone();
        let key = (organism.to_owned(), stage.to_owned(), identifier.to_owned());
        tokio::task::spawn_blocking(move || {
            let value = serde_json::to_string(&values).map_err(|e| e.to_string())?;
            this.connection
                .lock()
                .map_err(|e| e.to_string())?
                .execute(
                    "INSERT OR REPLACE INTO kegg_cache VALUES (?1, ?2, ?3, ?4, ?5)",
                    params![key.0, key.1, key.2, now(), value],
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
