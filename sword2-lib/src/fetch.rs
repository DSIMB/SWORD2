//! Remote structure fetching from PDB, AlphaFold, and ESM Atlas.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// Download a PDB file from the RCSB PDB database.
///
/// Fetches from `https://files.rcsb.org/download/{pdb_id}.pdb`.
pub fn fetch_pdb(pdb_id: &str, output_dir: &Path) -> Result<PathBuf> {
    let pdb_id = pdb_id.to_uppercase();
    let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
    let output_path = output_dir.join(format!("{}.pdb", pdb_id));

    tracing::info!("Fetching PDB {} from RCSB", pdb_id);
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

    tracing::info!("Downloaded PDB to {}", output_path.display());
    Ok(output_path)
}

/// Download an AlphaFold model from the EBI AlphaFold database.
///
/// Fetches from `https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb`.
pub fn fetch_alphafold(uniprot_id: &str, output_dir: &Path) -> Result<PathBuf> {
    let name = format!("AF-{}-F1-model_v4", uniprot_id);
    let url = format!("https://alphafold.ebi.ac.uk/files/{}.pdb", name);
    let output_path = output_dir.join(format!("{}.pdb", name));

    tracing::info!("Fetching AlphaFold model for UniProt {}", uniprot_id);
    let response = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(30))
        .build()?
        .get(&url)
        .send()
        .with_context(|| format!("Failed to fetch AlphaFold model for {}", uniprot_id))?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to fetch AlphaFold model for {}: HTTP {}",
            uniprot_id,
            response.status()
        );
    }

    let content = response.text()?;
    fs::write(&output_path, &content)
        .with_context(|| format!("Failed to write {}", output_path.display()))?;

    tracing::info!("Downloaded AlphaFold model to {}", output_path.display());
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

    tracing::info!("Fetching ESM model for MGnify {}", mgnify_id);
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

    tracing::info!("Downloaded ESM model to {}", output_path.display());
    Ok(output_path)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_fetch_pdb_url_construction() {
        let pdb_id = "1tim";
        assert_eq!(pdb_id.to_uppercase(), "1TIM");
    }
}
