use reqwest::{Client, Response};
use std::time::Duration;

#[derive(Clone)]
pub struct Upstream {
    pub client: Client,
    pub uniprot: String,
    pub quickgo: String,
}
impl Upstream {
    pub fn new(uniprot: String, quickgo: String) -> Result<Self, reqwest::Error> {
        Ok(Self {
            client: Client::builder()
                .connect_timeout(Duration::from_secs(5))
                .timeout(Duration::from_secs(30))
                .user_agent("protein-tools/0.1.0")
                .build()?,
            uniprot,
            quickgo,
        })
    }
    pub async fn get(&self, url: &str, params: &[(&str, String)]) -> Result<Response, String> {
        // No automatic retries: each request has a fixed timeout and traffic budget.
        let response = self
            .client
            .get(url)
            .query(params)
            .header("Accept", "application/json")
            .send()
            .await
            .map_err(|e| {
                tracing::warn!(event="upstream_failure", kind="network", error=%e);
                "upstream network request failed".to_string()
            })?;
        if response.status().is_client_error() || response.status().is_server_error() {
            tracing::warn!(event="upstream_failure", status=%response.status());
        }
        Ok(response)
    }
}
