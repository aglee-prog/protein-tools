use protein_tools::{
    api::{App, MAX_CONCURRENT_LOOKUPS, router},
    cache::Cache,
    upstream::Upstream,
};
use std::{env, sync::Arc};
use tokio::sync::Semaphore;
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    tracing_subscriber::fmt()
        .json()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "protein_tools=info".into()),
        )
        .init();
    let path = env::var("CACHE_DB").unwrap_or_else(|_| "/data/cache.sqlite".into());
    let ttl = env::var("CACHE_TTL_DAYS")
        .unwrap_or_else(|_| "30".into())
        .parse()?;
    let app = App {
        cache: Cache::open(&path, ttl)?,
        upstream: Upstream::new(
            "https://rest.uniprot.org".into(),
            "https://www.ebi.ac.uk/QuickGO/services".into(),
        )?,
        permits: Arc::new(Semaphore::new(MAX_CONCURRENT_LOOKUPS)),
    };
    let address = env::var("BIND_ADDR").unwrap_or_else(|_| "0.0.0.0:8080".into());
    let listener = tokio::net::TcpListener::bind(&address).await?;
    tracing::info!(event="server_startup", %address, cache_db=%path, cache_ttl_days=ttl);
    axum::serve(listener, router(app))
        .with_graceful_shutdown(shutdown())
        .await?;
    Ok(())
}
async fn shutdown() {
    let mut terminate = tokio::signal::unix::signal(tokio::signal::unix::SignalKind::terminate())
        .expect("SIGTERM handler");
    tokio::select! { _ = tokio::signal::ctrl_c() => {}, _ = terminate.recv() => {} }
}
