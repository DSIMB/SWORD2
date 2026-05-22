//! Amino acid three-letter to one-letter code conversion.
//!
//! Ported from the Python PDB.py `dico_AA` dictionary, which includes
//! the 20 standard amino acids plus many non-standard/modified residues.

use std::collections::HashMap;
use std::sync::LazyLock;

/// The 20 standard three-letter amino acid codes.
pub const STANDARD_AA3: &[&str] = &[
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO",
    "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
];

/// Mapping from three-letter amino acid codes (standard + non-standard) to
/// one-letter codes. Ported from the Python `dico_AA` dictionary.
static AA_MAP: LazyLock<HashMap<&'static str, char>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    // 20 standard amino acids
    m.insert("ALA", 'A');
    m.insert("CYS", 'C');
    m.insert("ASP", 'D');
    m.insert("GLU", 'E');
    m.insert("PHE", 'F');
    m.insert("GLY", 'G');
    m.insert("HIS", 'H');
    m.insert("ILE", 'I');
    m.insert("LYS", 'K');
    m.insert("LEU", 'L');
    m.insert("MET", 'M');
    m.insert("ASN", 'N');
    m.insert("PRO", 'P');
    m.insert("GLN", 'Q');
    m.insert("ARG", 'R');
    m.insert("SER", 'S');
    m.insert("THR", 'T');
    m.insert("VAL", 'V');
    m.insert("TRP", 'W');
    m.insert("TYR", 'Y');
    // Non-standard / modified residues
    m.insert("PAQ", 'Y');
    m.insert("AGM", 'R');
    m.insert("PR3", 'C');
    m.insert("DOH", 'D');
    m.insert("CCS", 'C');
    m.insert("GSC", 'G');
    m.insert("GHG", 'Q');
    m.insert("OAS", 'S');
    m.insert("MIS", 'S');
    m.insert("SIN", 'D');
    m.insert("TPL", 'W');
    m.insert("SAC", 'S');
    m.insert("4HT", 'W');
    m.insert("FGP", 'C');
    m.insert("HSO", 'H');
    m.insert("LYZ", 'K');
    m.insert("FGL", 'S');
    m.insert("PRS", 'P');
    m.insert("DCY", 'C');
    m.insert("LYM", 'K');
    m.insert("GPL", 'K');
    m.insert("PYX", 'C');
    m.insert("PCC", 'P');
    m.insert("EHP", 'F');
    m.insert("CHG", 'A');
    m.insert("TPO", 'T');
    m.insert("DAS", 'D');
    m.insert("AYA", 'A');
    m.insert("TYN", 'Y');
    m.insert("SVA", 'S');
    m.insert("SCY", 'C');
    m.insert("BNN", 'A');
    m.insert("5HP", 'E');
    m.insert("HAR", 'R');
    m.insert("IAS", 'D');
    m.insert("SNC", 'C');
    m.insert("AHB", 'N');
    m.insert("PTR", 'Y');
    m.insert("PHI", 'F');
    m.insert("NPH", 'C');
    m.insert("PHL", 'F');
    m.insert("SNN", 'D');
    m.insert("A66", 'A');
    m.insert("TYB", 'Y');
    m.insert("PHD", 'D');
    m.insert("MAA", 'A');
    m.insert("APN", 'A');
    m.insert("TYY", 'Y');
    m.insert("TYT", 'Y');
    m.insert("TIH", 'A');
    m.insert("TRG", 'K');
    m.insert("CXM", 'M');
    m.insert("DIV", 'V');
    m.insert("TYS", 'Y');
    m.insert("DTH", 'T');
    m.insert("MLE", 'L');
    m.insert("CME", 'C');
    m.insert("SHR", 'K');
    m.insert("OCY", 'C');
    m.insert("DTY", 'Y');
    m.insert("2AS", 'D');
    m.insert("AEI", 'T');
    m.insert("DTR", 'W');
    m.insert("OCS", 'C');
    m.insert("CMT", 'C');
    m.insert("BET", 'G');
    m.insert("NLP", 'L');
    m.insert("LLY", 'K');
    m.insert("SCH", 'C');
    m.insert("CEA", 'C');
    m.insert("LLP", 'K');
    m.insert("TRF", 'W');
    m.insert("HMR", 'R');
    m.insert("TYI", 'Y');
    m.insert("TRO", 'W');
    m.insert("NLE", 'L');
    m.insert("BMT", 'T');
    m.insert("BUC", 'C');
    m.insert("PEC", 'C');
    m.insert("BUG", 'L');
    m.insert("SCS", 'C');
    m.insert("NLN", 'L');
    m.insert("MHO", 'M');
    m.insert("CSO", 'C');
    m.insert("FTR", 'W');
    m.insert("DLE", 'L');
    m.insert("TRN", 'W');
    m.insert("CSE", 'C');
    m.insert("CSD", 'A');
    m.insert("OMT", 'M');
    m.insert("CSA", 'C');
    m.insert("DSP", 'D');
    m.insert("CSB", 'C');
    m.insert("DSN", 'S');
    m.insert("SHC", 'C');
    m.insert("CSX", 'C');
    m.insert("YCM", 'C');
    m.insert("CSZ", 'C');
    m.insert("TRQ", 'W');
    m.insert("CSW", 'C');
    m.insert("EFC", 'C');
    m.insert("CSP", 'C');
    m.insert("CSS", 'C');
    m.insert("CSR", 'C');
    m.insert("CZZ", 'C');
    m.insert("MSO", 'M');
    m.insert("BTR", 'W');
    m.insert("HLU", 'L');
    m.insert("MGN", 'Q');
    m.insert("HTI", 'C');
    m.insert("TYQ", 'Y');
    m.insert("4IN", 'W');
    m.insert("M3L", 'K');
    m.insert("C5C", 'C');
    m.insert("HTR", 'W');
    m.insert("MPQ", 'G');
    m.insert("KCX", 'K');
    m.insert("GLH", 'E');
    m.insert("DIL", 'I');
    m.insert("ACA", 'A');
    m.insert("NEM", 'H');
    m.insert("5CS", 'C');
    m.insert("LYX", 'K');
    m.insert("DVA", 'V');
    m.insert("ACL", 'R');
    m.insert("GLX", 'Z');
    m.insert("MLZ", 'K');
    m.insert("GLZ", 'G');
    m.insert("SME", 'M');
    m.insert("SMC", 'C');
    m.insert("DLY", 'K');
    m.insert("NEP", 'H');
    m.insert("BCS", 'C');
    m.insert("ASQ", 'D');
    m.insert("SET", 'S');
    m.insert("SEP", 'S');
    m.insert("ASX", 'B');
    m.insert("DGN", 'Q');
    m.insert("DGL", 'E');
    m.insert("MHS", 'H');
    m.insert("SEG", 'A');
    m.insert("ASB", 'D');
    m.insert("ASA", 'D');
    m.insert("SEC", 'C');
    m.insert("SEB", 'S');
    m.insert("ASK", 'D');
    m.insert("GGL", 'E');
    m.insert("ASI", 'N');
    m.insert("SEL", 'S');
    m.insert("CGU", 'E');
    m.insert("C6C", 'C');
    m.insert("ASL", 'D');
    m.insert("LTR", 'W');
    m.insert("CLD", 'S');
    m.insert("CLE", 'L');
    m.insert("GMA", 'E');
    m.insert("1LU", 'L');
    m.insert("CLB", 'S');
    m.insert("MVA", 'V');
    m.insert("S1H", 'S');
    m.insert("DNP", 'A');
    m.insert("SAR", 'G');
    m.insert("FME", 'M');
    m.insert("ALO", 'T');
    m.insert("ALM", 'A');
    m.insert("LEF", 'L');
    m.insert("MEN", 'N');
    m.insert("TPQ", 'Y');
    m.insert("NMC", 'G');
    m.insert("SBD", 'S');
    m.insert("ALY", 'K');
    m.insert("MME", 'M');
    m.insert("GL3", 'G');
    m.insert("ALS", 'C');
    m.insert("SBL", 'S');
    m.insert("2MR", 'R');
    m.insert("CAY", 'C');
    m.insert("3AH", 'H');
    m.insert("DPR", 'P');
    m.insert("CAS", 'C');
    m.insert("NC1", 'S');
    m.insert("HYP", 'P');
    m.insert("FLA", 'A');
    m.insert("LCX", 'K');
    m.insert("MSE", 'M');
    m.insert("IYR", 'Y');
    m.insert("DPN", 'F');
    m.insert("BAL", 'A');
    m.insert("CAF", 'C');
    m.insert("MSA", 'G');
    m.insert("AIB", 'A');
    m.insert("HIP", 'H');
    m.insert("CYQ", 'C');
    m.insert("PCA", 'E');
    m.insert("DAL", 'A');
    m.insert("BFD", 'D');
    m.insert("DAH", 'F');
    m.insert("HIC", 'H');
    m.insert("CYG", 'C');
    m.insert("DAR", 'R');
    m.insert("CYD", 'C');
    m.insert("IIL", 'I');
    m.insert("CYM", 'C');
    m.insert("CYL", 'C');
    m.insert("CY3", 'C');
    m.insert("CY1", 'C');
    m.insert("HAC", 'A');
    m.insert("143", 'C');
    m.insert("DHI", 'H');
    m.insert("CY4", 'C');
    m.insert("YOF", 'Y');
    m.insert("HPQ", 'F');
    m.insert("SOC", 'C');
    m.insert("DHA", 'A');
    m.insert("2LU", 'L');
    m.insert("MLY", 'K');
    m.insert("TRW", 'W');
    m.insert("STY", 'Y');
    m.insert("MCL", 'K');
    m.insert("BHD", 'D');
    m.insert("NRQ", 'Y');
    m.insert("ARM", 'R');
    m.insert("PRR", 'A');
    m.insert("ARO", 'R');
    // Protonated histidines (CHARMM)
    m.insert("HSE", 'H');
    m.insert("HSP", 'H');
    m.insert("HSD", 'H');
    // Non-standard from AA3 list
    m.insert("5HP", 'E');
    m.insert("ABA", 'A');
    m.insert("RON", 'X');
    m.insert("3GA", 'X');
    m
});

/// Convert a three-letter amino acid code to a one-letter code.
///
/// Returns `Some(char)` for known residues, `None` for unknown ones.
pub fn three_to_one(three: &str) -> Option<char> {
    AA_MAP.get(three.trim()).copied()
}

/// Returns `true` if the given three-letter code is one of the 20 standard amino acids.
pub fn is_standard(three: &str) -> bool {
    STANDARD_AA3.contains(&three.trim())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_amino_acids() {
        assert_eq!(three_to_one("ALA"), Some('A'));
        assert_eq!(three_to_one("TRP"), Some('W'));
        assert_eq!(three_to_one("TYR"), Some('Y'));
    }

    #[test]
    fn test_non_standard() {
        assert_eq!(three_to_one("MSE"), Some('M'));
        assert_eq!(three_to_one("HSE"), Some('H'));
    }

    #[test]
    fn test_unknown() {
        assert_eq!(three_to_one("ZZZ"), None);
    }

    #[test]
    fn test_is_standard() {
        assert!(is_standard("ALA"));
        assert!(!is_standard("MSE"));
    }
}
