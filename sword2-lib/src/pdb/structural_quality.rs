use anyhow::Result;
use serde::Serialize;

use super::{amino_acids, Atom, Chain, Residue};

pub const STRUCTURAL_ELIGIBILITY_POLICY: &str = "strict_complete_backbone_v1";
pub const REQUIRED_BACKBONE_ATOMS: [&str; 4] = ["N", "CA", "C", "O"];

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct IncompleteBackboneResidue {
    pub author_residue_number: i32,
    pub chain_id: String,
    pub missing_atoms: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct StructuralQualityReport {
    pub candidate_residue_count: usize,
    pub chain_id: String,
    pub complete_backbone_residue_count: usize,
    pub eligible: bool,
    pub incomplete_residues: Vec<IncompleteBackboneResidue>,
    pub policy: &'static str,
    pub reason_code: Option<&'static str>,
    pub schema_version: u32,
    pub structural_coverage: f64,
}

pub fn accepts_backbone_alt_loc(alt_loc: char) -> bool {
    matches!(alt_loc, ' ' | 'A')
}

pub fn is_sword_candidate_residue(residue: &Residue) -> bool {
    residue.atoms.iter().any(|atom| !atom.is_hetatm)
        && amino_acids::is_standard(&residue.name)
        && residue.icode == ' '
}

fn accepted_required_atom(atom: &Atom, name: &str) -> bool {
    !atom.is_hetatm
        && accepts_backbone_alt_loc(atom.alt_loc)
        && atom.name.trim() == name
        && atom.coord.x.is_finite()
        && atom.coord.y.is_finite()
        && atom.coord.z.is_finite()
}

pub fn missing_backbone_atoms(residue: &Residue) -> Vec<String> {
    REQUIRED_BACKBONE_ATOMS
        .iter()
        .filter(|name| {
            !residue
                .atoms
                .iter()
                .any(|atom| accepted_required_atom(atom, name))
        })
        .map(|name| (*name).to_owned())
        .collect()
}

pub fn inspect_structural_quality(chain: &Chain) -> StructuralQualityReport {
    let candidates: Vec<&Residue> = chain
        .residues
        .iter()
        .filter(|residue| is_sword_candidate_residue(residue))
        .collect();
    let incomplete_residues: Vec<_> = candidates
        .iter()
        .filter_map(|residue| {
            let missing_atoms = missing_backbone_atoms(residue);
            (!missing_atoms.is_empty()).then(|| IncompleteBackboneResidue {
                author_residue_number: residue.seq_num,
                chain_id: residue.chain_id.to_string(),
                missing_atoms,
            })
        })
        .collect();
    let complete_backbone_residue_count = candidates.len() - incomplete_residues.len();

    StructuralQualityReport {
        candidate_residue_count: candidates.len(),
        chain_id: chain.id.to_string(),
        complete_backbone_residue_count,
        eligible: !candidates.is_empty() && incomplete_residues.is_empty(),
        incomplete_residues,
        policy: STRUCTURAL_ELIGIBILITY_POLICY,
        reason_code: candidates
            .is_empty()
            .then_some("empty_candidate_population"),
        schema_version: 1,
        structural_coverage: if candidates.is_empty() {
            0.0
        } else {
            complete_backbone_residue_count as f64 / candidates.len() as f64
        },
    }
}

impl StructuralQualityReport {
    pub fn canonical_json_bytes(&self) -> Result<Vec<u8>> {
        let mut bytes = serde_json::to_vec(self)?;
        bytes.push(b'\n');
        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pdb::{Atom, Chain, Point3D, Residue};

    fn atom(name: &str, alt_loc: char, xyz: [f64; 3], is_hetatm: bool) -> Atom {
        Atom::new(
            1,
            name,
            alt_loc,
            "ALA",
            'A',
            42,
            ' ',
            Point3D::new(xyz[0], xyz[1], xyz[2]),
            1.0,
            10.0,
            "C",
            "",
            is_hetatm,
        )
    }

    fn complete_residue() -> Residue {
        let mut residue = Residue::new("ALA", 42, ' ', 'A');
        for name in REQUIRED_BACKBONE_ATOMS {
            residue.atoms.push(atom(name, ' ', [1.0, 2.0, 3.0], false));
        }
        residue
    }

    #[test]
    fn complete_standard_residue_is_eligible() {
        let mut chain = Chain::new('A');
        chain.residues.push(complete_residue());
        let report = inspect_structural_quality(&chain);
        assert!(report.eligible);
        assert_eq!(report.candidate_residue_count, 1);
        assert_eq!(report.complete_backbone_residue_count, 1);
        assert_eq!(report.structural_coverage, 1.0);
        assert!(report.incomplete_residues.is_empty());
    }

    #[test]
    fn missing_and_nonfinite_atoms_are_reported_in_canonical_order() {
        let mut residue = complete_residue();
        residue
            .atoms
            .retain(|atom| !matches!(atom.name.trim(), "CA" | "O"));
        residue
            .atoms
            .push(atom("O", ' ', [f64::NAN, 0.0, 0.0], false));
        let mut chain = Chain::new('A');
        chain.residues.push(residue);
        let report = inspect_structural_quality(&chain);
        assert!(!report.eligible);
        assert_eq!(
            report.incomplete_residues[0].missing_atoms,
            vec!["CA".to_string(), "O".to_string()],
        );
    }

    #[test]
    fn only_blank_or_a_altloc_can_supply_required_atoms() {
        let mut residue = complete_residue();
        residue.atoms.retain(|atom| atom.name.trim() != "CA");
        residue.atoms.push(atom("CA", 'B', [1.0, 2.0, 3.0], false));
        let mut chain = Chain::new('A');
        chain.residues.push(residue);
        assert!(!inspect_structural_quality(&chain).eligible);
        chain.residues[0]
            .atoms
            .push(atom("CA", 'A', [1.0, 2.0, 3.0], false));
        assert!(inspect_structural_quality(&chain).eligible);
    }

    #[test]
    fn residues_outside_sword_cleaning_do_not_enter_coverage() {
        let mut chain = Chain::new('A');
        let mut insertion = complete_residue();
        insertion.icode = 'A';
        let mut nonstandard = complete_residue();
        nonstandard.name = "MSE".into();
        let mut hetatm = complete_residue();
        for atom in &mut hetatm.atoms {
            atom.is_hetatm = true;
        }
        chain.residues.extend([insertion, nonstandard, hetatm]);
        let report = inspect_structural_quality(&chain);
        assert_eq!(report.candidate_residue_count, 0);
        assert_eq!(report.structural_coverage, 0.0);
        assert_eq!(
            report.reason_code.as_deref(),
            Some("empty_candidate_population")
        );
    }

    #[test]
    fn report_bytes_are_canonical_and_newline_terminated() {
        let mut chain = Chain::new('A');
        chain.residues.push(complete_residue());
        let bytes = inspect_structural_quality(&chain)
            .canonical_json_bytes()
            .unwrap();
        assert_eq!(bytes.last(), Some(&b'\n'));
        assert_eq!(
            bytes,
            b"{\"candidate_residue_count\":1,\"chain_id\":\"A\",\"complete_backbone_residue_count\":1,\"eligible\":true,\"incomplete_residues\":[],\"policy\":\"strict_complete_backbone_v1\",\"reason_code\":null,\"schema_version\":1,\"structural_coverage\":1.0}\n",
        );
    }
}
