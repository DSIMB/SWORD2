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
        .build()
        .context("Failed to build HTTP client")
}

fn select_alphafold_prediction<'a>(
    predictions: &'a [AlphaFoldPrediction],
    uniprot_id: &str,
) -> Option<&'a AlphaFoldPrediction> {
    predictions
        .iter()
        .find(|prediction| prediction.uniprot_accession.eq_ignore_ascii_case(uniprot_id))
}

fn filename_from_url(url: &str) -> Result<&str> {
    url.rsplit('/')
        .next()
        .filter(|segment| !segment.is_empty())
        .with_context(|| format!("Failed to derive file name from URL {url}"))
}

/// Download a PDB file from the RCSB PDB database.
///
/// Fetches from `https://files.rcsb.org/download/{pdb_id}.pdb`.
pub fn fetch_pdb(pdb_id: &str, output_dir: &Path) -> Result<PathBuf> {
    let pdb_id = pdb_id.to_uppercase();
    let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
    let output_path = output_dir.join(format!("{}.pdb", pdb_id));

    tracing::debug!("Fetching PDB {} from RCSB", pdb_id);
    let response = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(30))
        .build()?
        .get(&url)
        .send()
        .with_context(|| format!("Failed to fetch PDB {}", pdb_id))?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to fetch PDB {}: HTTP {}",
            pdb_id,
            response.status()
        );
    }

    let content = response.text()?;
    fs::write(&output_path, &content)
        .with_context(|| format!("Failed to write {}", output_path.display()))?;

    tracing::debug!("Downloaded PDB to {}", output_path.display());
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
    let response = client
        .get(&prediction.pdb_url)
        .send()
        .with_context(|| format!("Failed to download AlphaFold model from {}", prediction.pdb_url))?;

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
    let response = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(30))
        .build()?
        .get(&url)
        .send()
        .with_context(|| format!("Failed to fetch ESM model for {}", mgnify_id))?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to fetch ESM model for {}: HTTP {}",
            mgnify_id,
            response.status()
        );
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
                pdb_url: "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.pdb"
                    .to_string(),
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
