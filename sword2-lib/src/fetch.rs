//! Remote structure fetching from PDB, AlphaFold, and ESM Atlas.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct AlphaFoldPrediction {
    #[serde(rename = "uniprotAccession")]
    uniprot_accession: String,
    #[serde(rename = "pdbUrl")]
    pdb_url: String,
}

fn build_http_client() -> Result<reqwest::blocking::Client> {
    reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(30))
        // A User-Agent is required by some upstreams: the EBI AlphaFold API
        // returns HTTP 403 for requests without one. reqwest sends none by
        // default, so set an explicit identifying UA for all fetchers.
        .user_agent(concat!("SWORD2/", env!("CARGO_PKG_VERSION")))
        .build()
        .context("Failed to build HTTP client")
}

fn select_alphafold_prediction<'a>(
    predictions: &'a [AlphaFoldPrediction],
    uniprot_id: &str,
) -> Option<&'a AlphaFoldPrediction> {
    predictions.iter().find(|prediction| {
        prediction
            .uniprot_accession
            .eq_ignore_ascii_case(uniprot_id)
    })
}

fn filename_from_url(url: &str) -> Result<&str> {
    url.rsplit('/')
        .next()
        .filter(|segment| !segment.is_empty())
        .with_context(|| format!("Failed to derive file name from URL {url}"))
}

/// Download a structure from the RCSB PDB database.
///
/// By default fetches mmCIF (`.cif`), which is the current RCSB standard and
/// works for all entries including large structures. Pass `prefer_pdb = true`
/// to request the legacy `.pdb` format instead; if that 404s, the error is
/// surfaced directly (no silent fallback).
pub fn fetch_pdb(pdb_id: &str, output_dir: &Path, prefer_pdb: bool) -> Result<PathBuf> {
    let pdb_id = pdb_id.to_uppercase();
    let client = build_http_client()?;

    if prefer_pdb {
        let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
        tracing::debug!("Fetching {} in legacy PDB format from RCSB", pdb_id);
        let response = client
            .get(&url)
            .send()
            .with_context(|| format!("Failed to fetch PDB {}", pdb_id))?;
        if !response.status().is_success() {
            anyhow::bail!(
                "Failed to fetch {} in PDB format: HTTP {} \
                 (try without --legacy-pdb to use mmCIF instead)",
                pdb_id,
                response.status()
            );
        }
        let output_path = output_dir.join(format!("{}.pdb", pdb_id));
        fs::write(&output_path, response.text()?)
            .with_context(|| format!("Failed to write {}", output_path.display()))?;
        tracing::debug!("Downloaded PDB to {}", output_path.display());
        return Ok(output_path);
    }

    // Default: mmCIF
    let url = format!("https://files.rcsb.org/download/{}.cif", pdb_id);
    tracing::debug!("Fetching {} in mmCIF format from RCSB", pdb_id);
    let response = client
        .get(&url)
        .send()
        .with_context(|| format!("Failed to fetch mmCIF for {}", pdb_id))?;
    if !response.status().is_success() {
        anyhow::bail!("Failed to fetch {}: HTTP {}", pdb_id, response.status());
    }
    let output_path = output_dir.join(format!("{}.cif", pdb_id));
    fs::write(&output_path, response.text()?)
        .with_context(|| format!("Failed to write {}", output_path.display()))?;
    tracing::debug!("Downloaded mmCIF to {}", output_path.display());
    Ok(output_path)
}

/// Download an AlphaFold model from the EBI AlphaFold database.
///
/// Resolves the latest model URL from `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`.
pub fn fetch_alphafold(uniprot_id: &str, output_dir: &Path) -> Result<PathBuf> {
    let uniprot_id = uniprot_id.trim().to_ascii_uppercase();
    let metadata_url = format!("https://alphafold.ebi.ac.uk/api/prediction/{}", uniprot_id);
    let client = build_http_client()?;

    tracing::debug!("Fetching AlphaFold model for UniProt {}", uniprot_id);
    let response = client
        .get(&metadata_url)
        .header(reqwest::header::ACCEPT, "application/json")
        .send()
        .with_context(|| format!("Failed to fetch AlphaFold model for {}", uniprot_id))?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to fetch AlphaFold metadata for {}: HTTP {}",
            uniprot_id,
            response.status()
        );
    }

    let metadata = response.text()?;
    let predictions: Vec<AlphaFoldPrediction> = serde_json::from_str(&metadata)
        .with_context(|| format!("Failed to parse AlphaFold API response for {}", uniprot_id))?;
    let prediction = select_alphafold_prediction(&predictions, &uniprot_id).with_context(|| {
        format!(
            "AlphaFold API returned no prediction matching UniProt accession {}",
            uniprot_id
        )
    })?;

    let output_name = filename_from_url(&prediction.pdb_url)?;
    let output_path = output_dir.join(output_name);
    let response = client.get(&prediction.pdb_url).send().with_context(|| {
        format!(
            "Failed to download AlphaFold model from {}",
            prediction.pdb_url
        )
    })?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to download AlphaFold model for {}: HTTP {}",
            uniprot_id,
            response.status()
        );
    }

    let content = response.text()?;
    fs::write(&output_path, &content)
        .with_context(|| format!("Failed to write {}", output_path.display()))?;

    tracing::debug!("Downloaded AlphaFold model to {}", output_path.display());
    Ok(output_path)
}

/// Download an ESM Metagenomic Atlas model.
///
/// Fetches from `https://api.esmatlas.com/fetchPredictedStructure/{mgnify_id}`.
pub fn fetch_esm(mgnify_id: &str, output_dir: &Path) -> Result<PathBuf> {
    let url = format!(
        "https://api.esmatlas.com/fetchPredictedStructure/{}",
        mgnify_id
    );
    let output_path = output_dir.join(format!("{}.pdb", mgnify_id));

    tracing::debug!("Fetching ESM model for MGnify {}", mgnify_id);
    let response = build_http_client()?
        .get(&url)
        .send()
        .with_context(|| format!("Failed to fetch ESM model for {}", mgnify_id))?;

    if !response.status().is_success() {
        let status = response.status();
        // ESM Atlas returns 403 (not 404) for unknown IDs. The most common
        // mistake is passing a MGnify *study* accession (MGYS…); the structure
        // API only indexes *protein* accessions (MGYP…).
        if status == reqwest::StatusCode::FORBIDDEN || status == reqwest::StatusCode::NOT_FOUND {
            anyhow::bail!(
                "Failed to fetch ESM model for {}: HTTP {}. \
                 Ensure this is a valid MGnify protein accession (MGYP…); \
                 study accessions (MGYS…) are not valid ESM Atlas structure IDs.",
                mgnify_id,
                status
            );
        }
        anyhow::bail!("Failed to fetch ESM model for {}: HTTP {}", mgnify_id, status);
    }

    let content = response.text()?;
    fs::write(&output_path, &content)
        .with_context(|| format!("Failed to write {}", output_path.display()))?;

    tracing::debug!("Downloaded ESM model to {}", output_path.display());
    Ok(output_path)
}

#[cfg(test)]
mod tests {
    use super::{filename_from_url, select_alphafold_prediction, AlphaFoldPrediction};

    #[test]
    fn test_fetch_pdb_url_construction() {
        let pdb_id = "1tim";
        assert_eq!(pdb_id.to_uppercase(), "1TIM");
    }

    #[test]
    fn test_select_alphafold_prediction_matches_requested_accession() {
        let predictions = vec![
            AlphaFoldPrediction {
                uniprot_accession: "Q5VSL9".to_string(),
                pdb_url: "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.pdb".to_string(),
            },
            AlphaFoldPrediction {
                uniprot_accession: "Q5VSL9-4".to_string(),
                pdb_url: "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-4-F1-model_v6.pdb"
                    .to_string(),
            },
        ];

        let prediction = select_alphafold_prediction(&predictions, "q5vsl9-4").unwrap();
        assert_eq!(prediction.uniprot_accession, "Q5VSL9-4");
    }

    #[test]
    fn test_filename_from_url_uses_latest_version_name() {
        let url = "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.pdb";
        assert_eq!(filename_from_url(url).unwrap(), "AF-Q5VSL9-F1-model_v6.pdb");
    }
}
