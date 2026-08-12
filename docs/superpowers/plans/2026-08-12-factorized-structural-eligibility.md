# Factorized Structural Eligibility Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add strict complete-backbone eligibility and explicit scientific abstention to the frozen factorized ranker, then evaluate it only on a predeclared input-eligible subset while retaining all 663 CATH identities in the audit denominator.

**Architecture:** A new Rust structural-quality module owns the single residue/atom contract used by normal runs, the read-only inspector, DSSP alternate-location selection, and factorized preflight. A canonical Python eligibility artifact freezes input-only decisions before any replacement evaluation; locked collection and validation then accept only those exact abstentions and compute conditional comparisons on the same eligible identities. Existing model weights and ordinary legacy selection remain unchanged.

**Tech Stack:** Rust 2021, serde/serde_json, clap, Python 3.12, pytest, pandas, deterministic canonical JSON, SHA-256, Git-bound runtime freezes.

**Execution mode:** Inline execution in the current workspace, as requested. Do not create a worktree or dispatch subagents. Preserve unrelated dirty/untracked files. Parallelize only independent verification and accuracy processes; paired resource measurements remain sequential.

---

## Scientific and operational guardrails

- Do not change `benchmark/models/factorized_count_v1.json`, `benchmark/models/factorized_candidate_v1.json`, generated Rust tree bytes, feature names, or tree thresholds.
- Do not filter, impute, or renumber incomplete residues for factorized model version 1.
- Do not inspect another CATH-663 score or competitor output until the replacement runtime, eligibility-inspection intent, and input-only eligibility artifact are committed.
- Keep the two previous failed ledgers because they are audit evidence, not obsolete scratch data.
- Delete temporary development output after its hashes/counts are recorded.
- Use CATH-17287 only for development validation.
- An incomplete chain may still return the legacy partition to users, but its status must state that no factorized decision occurred.
- One missing required atom is sufficient for abstention. Version 1 has no learned or percentage threshold.
- The canonical required-atom order is `N`, `CA`, `C`, `O`.
- Accuracy runs use 32 concurrent SWORD processes with one Rust thread each. Resource runs use exactly one process at a time in the frozen counterbalanced order.

## File map

### New files

- `sword2-lib/src/pdb/structural_quality.rs` — authoritative residue eligibility, report schema, and canonical serialization.
- `benchmark/factorized_ranker/eligibility.py` — canonical report/manifest validation and deterministic parallel input inspection.
- `benchmark/freeze_factorized_eligibility.py` — `create`/`verify` CLI for the input-only eligibility freeze.
- `benchmark/tests/test_factorized_eligibility.py` — eligibility report, artifact, tamper, and identity-prefix tests.
- Generated only after source freeze: `benchmark/models/factorized_ranker_v1_runtime_v2.json`.
- Generated only after a committed eligibility intent: `benchmark/models/factorized_ranker_v1_structural_eligibility.json`.
- Generated only before replacement collection: attempt-3 eligibility/evaluation intents and result paths; names are frozen in Task 9.

### Modified Rust files

- `sword2-lib/src/pdb/mod.rs` — export structural-quality interfaces.
- `sword2-lib/src/pdb/writer.rs` — reuse the authoritative base candidate-residue predicate.
- `sword2-lib/src/dssp/backbone.rs` — reuse the authoritative alternate-location predicate.
- `sword2-lib/src/sword/factorized_ranker/mod.rs` — stable abstention error/status code.
- `sword2-lib/src/sword/mod.rs` — block learned inference before feature/model execution while retaining legacy output.
- `sword2-cli/src/main.rs` — read-only inspection subcommand, shared model/chain selection, and normal-run quality report.

### Modified Python files

- `benchmark/factorized_ranker/runtime_freeze.py` — runtime schema v2 and new evidence-tool hashes while retaining schema-v1 parsing.
- `benchmark/run_benchmark.py` — eligibility-bound locked preflight, expected abstentions, conditional resource set, and 32-worker accuracy collection.
- `benchmark/evaluate_factorized_acceptance.py` — eligibility-aware coverage, conditional gates, full-denominator coverage reporting, and fail-closed promotion.
- `benchmark/tests/test_factorized_runtime_freeze.py` — schema-v2 evidence graph.
- `benchmark/tests/test_run_benchmark.py` — expected-abstention and bounded-parallel collection behavior.
- `benchmark/tests/test_factorized_acceptance.py` — revised coverage/evaluation tamper matrix.

## Task 1: Build the authoritative Rust structural-quality contract

**Files:**
- Create: `sword2-lib/src/pdb/structural_quality.rs`
- Modify: `sword2-lib/src/pdb/mod.rs`
- Modify: `sword2-lib/src/pdb/writer.rs:1-67`
- Modify: `sword2-lib/src/dssp/backbone.rs:35-110`

- [ ] **Step 1: Write the failing structural-quality tests**

Create the new module with tests that construct `Residue`/`Atom` fixtures and call the initially absent production API:

```rust
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
        residue.atoms.retain(|atom| !matches!(atom.name.trim(), "CA" | "O"));
        residue.atoms.push(atom("O", ' ', [f64::NAN, 0.0, 0.0], false));
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
        chain.residues[0].atoms.push(atom("CA", 'A', [1.0, 2.0, 3.0], false));
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
        for atom in &mut hetatm.atoms { atom.is_hetatm = true; }
        chain.residues.extend([insertion, nonstandard, hetatm]);
        let report = inspect_structural_quality(&chain);
        assert_eq!(report.candidate_residue_count, 0);
        assert_eq!(report.structural_coverage, 0.0);
        assert_eq!(report.reason_code.as_deref(), Some("empty_candidate_population"));
    }

    #[test]
    fn report_bytes_are_canonical_and_newline_terminated() {
        let mut chain = Chain::new('A');
        chain.residues.push(complete_residue());
        let bytes = inspect_structural_quality(&chain).canonical_json_bytes().unwrap();
        assert_eq!(bytes.last(), Some(&b'\n'));
        assert_eq!(
            bytes,
            b"{\"candidate_residue_count\":1,\"chain_id\":\"A\",\"complete_backbone_residue_count\":1,\"eligible\":true,\"incomplete_residues\":[],\"policy\":\"strict_complete_backbone_v1\",\"reason_code\":null,\"schema_version\":1,\"structural_coverage\":1.0}\n",
        );
    }
}
```

- [ ] **Step 2: Run the focused test and capture RED**

Run:

```bash
cargo test -p sword2-lib pdb::structural_quality -- --nocapture
```

Expected: compilation fails because `REQUIRED_BACKBONE_ATOMS`, `inspect_structural_quality`, and report types do not exist.

- [ ] **Step 3: Implement the minimal typed contract**

Implement these public interfaces in `structural_quality.rs`:

```rust
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
        .filter(|name| !residue.atoms.iter().any(|atom| accepted_required_atom(atom, name)))
        .map(|name| (*name).to_owned())
        .collect()
}

pub fn inspect_structural_quality(chain: &Chain) -> StructuralQualityReport {
    let candidates: Vec<&Residue> = chain.residues.iter().filter(|r| is_sword_candidate_residue(r)).collect();
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
    let complete = candidates.len() - incomplete_residues.len();
    StructuralQualityReport {
        candidate_residue_count: candidates.len(),
        chain_id: chain.id.to_string(),
        complete_backbone_residue_count: complete,
        eligible: !candidates.is_empty() && incomplete_residues.is_empty(),
        incomplete_residues,
        policy: STRUCTURAL_ELIGIBILITY_POLICY,
        reason_code: candidates.is_empty().then_some("empty_candidate_population"),
        schema_version: 1,
        structural_coverage: if candidates.is_empty() { 0.0 } else { complete as f64 / candidates.len() as f64 },
    }
}

impl StructuralQualityReport {
    pub fn canonical_json_bytes(&self) -> Result<Vec<u8>> {
        let mut bytes = serde_json::to_vec(self)?;
        bytes.push(b'\n');
        Ok(bytes)
    }
}
```

Keep fields in lexical order so serde output is canonical without an untyped intermediate map.

- [ ] **Step 4: Make cleaning and DSSP share the predicates**

Export the module in `pdb/mod.rs` and replace the three duplicated cleaning conditions in both loops of `clean_chain_for_sword` with:

```rust
if !super::structural_quality::is_sword_candidate_residue(residue) {
    continue;
}
```

In `dssp/backbone.rs`, replace:

```rust
if altloc != ' ' && altloc != 'A' {
    continue;
}
```

with:

```rust
if !crate::pdb::structural_quality::accepts_backbone_alt_loc(altloc) {
    continue;
}
```

- [ ] **Step 5: Run focused and neighboring tests GREEN**

Run:

```bash
cargo test -p sword2-lib pdb::structural_quality -- --nocapture
cargo test -p sword2-lib pdb:: -- --nocapture
cargo test -p sword2-lib dssp:: -- --nocapture
```

Expected: all pass; no legacy cleaning fixture changes.

- [ ] **Step 6: Commit Task 1**

```bash
git add sword2-lib/src/pdb/structural_quality.rs \
  sword2-lib/src/pdb/mod.rs \
  sword2-lib/src/pdb/writer.rs \
  sword2-lib/src/dssp/backbone.rs
git diff --cached --check
git commit -m "feat: inspect factorized structural eligibility"
```

## Task 2: Add the read-only inspector and normal-run quality artifact

**Files:**
- Modify: `sword2-cli/src/main.rs:1-180,523-700,930-950,1100-1185`

- [ ] **Step 1: Write failing CLI parsing and byte-identity tests**

Add tests that require the new subcommand and shared selection helper:

```rust
#[test]
fn parses_factorized_eligibility_inspector() {
    let cli = Cli::try_parse_from([
        "sword2", "inspect-factorized-eligibility",
        "--input", "fixture.pdb", "--chain", "B", "--nmr-model", "2",
    ]).unwrap();
    assert!(matches!(
        cli.command,
        Some(Commands::InspectFactorizedEligibility(InspectFactorizedEligibilityArgs {
            chain: Some('B'), nmr_model: 2, ..
        }))
    ));
}

#[test]
fn inspector_and_normal_report_use_identical_bytes() {
    let chain = complete_chain_fixture();
    let direct = pdb::structural_quality::inspect_structural_quality(&chain)
        .canonical_json_bytes().unwrap();
    assert_eq!(structural_quality_bytes(&chain).unwrap(), direct);
}
```

Also add a selection test proving model serial `2` and explicit chain `B` are honored and absent model/chain IDs fail.

- [ ] **Step 2: Run CLI tests and capture RED**

Run:

```bash
cargo test -p sword parses_factorized_eligibility_inspector -- --nocapture
```

Expected: compilation fails for the absent command, args, helper, and byte function.

- [ ] **Step 3: Implement the CLI interface and shared selection**

Add:

```rust
#[derive(Args, Debug)]
struct InspectFactorizedEligibilityArgs {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    chain: Option<char>,
    #[arg(long, default_value = "1")]
    nmr_model: i32,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Score(ScoreArgs),
    InspectFactorizedEligibility(InspectFactorizedEligibilityArgs),
}
```

Extract model/chain selection into one helper used by `process_entry` and the inspector:

```rust
fn select_chain<'a>(
    structure: &'a pdb::Structure,
    model_serial: i32,
    requested_chain: Option<char>,
) -> Result<&'a pdb::Chain> {
    let model = structure.get_model(model_serial)
        .ok_or_else(|| anyhow::anyhow!("Model {model_serial} not found"))?;
    if let Some(chain_id) = requested_chain {
        return model.get_chain(chain_id).ok_or_else(|| {
            let available = model.chain_ids().into_iter().map(|c| c.to_string()).collect::<Vec<_>>();
            anyhow::anyhow!("Chain '{chain_id}' not found. Available chains: {}", available.join(", "))
        });
    }
    model.chains.iter()
        .find(|chain| chain.residues.iter().any(pdb::structural_quality::is_sword_candidate_residue))
        .or_else(|| model.chains.first())
        .ok_or_else(|| anyhow::anyhow!("No chains found in model {model_serial}"))
}
```

The inspector parses once, selects once, computes once, writes only canonical bytes to stdout, and returns zero for either `eligible=true` or `eligible=false`:

```rust
fn run_inspect_factorized_eligibility(args: &InspectFactorizedEligibilityArgs) -> Result<()> {
    let structure = pdb::parse_pdb(&args.input)
        .with_context(|| format!("Failed to parse {}", args.input.display()))?;
    let chain = select_chain(&structure, args.nmr_model, args.chain)?;
    use std::io::Write as _;
    std::io::stdout().write_all(&structural_quality_bytes(chain)?)?;
    Ok(())
}
```

- [ ] **Step 4: Write `residue_quality.json` before normal pipeline execution**

After optional pLDDT filtering and before cleaning, compute the report and atomically write its exact bytes to `results_dir/residue_quality.json`. Use a sibling temporary file, `write_all`, `flush`, `sync_all`, and `rename`; delete only the exact temporary file after an error. Then call the unchanged cleaning function and legacy pipeline.

Do not use the report to filter the chain. Do not change a legacy selection or domain mapping.

- [ ] **Step 5: Run CLI and library regression tests GREEN**

```bash
cargo test -p sword -- --nocapture
cargo test -p sword2-lib pdb:: -- --nocapture
```

Expected: inspector tests pass, default CLI behavior remains factorized-off, and normal report bytes match inspector bytes.

- [ ] **Step 6: Commit Task 2**

```bash
git add sword2-cli/src/main.rs
git diff --cached --check
git commit -m "feat: expose structural eligibility inspection"
```

## Task 3: Enforce abstention before factorized feature/model execution

**Files:**
- Modify: `sword2-lib/src/sword/factorized_ranker/mod.rs:30-145,1680-1785`
- Modify: `sword2-lib/src/sword/mod.rs:110-190,245-410,750-835,1290-1360`

- [ ] **Step 1: Write failing status and transaction tests**

Add exact tests:

```rust
#[test]
fn structural_abstention_has_stable_error_code_and_status_bytes() {
    let error = FactorizedError::StructuralQualityAbstention;
    assert_eq!(error.code(), "structural_quality_abstention");
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("status.json");
    write_selector_status(&path, &SelectorStatus::factorized_fallback(&error, 0)).unwrap();
    assert_eq!(
        std::fs::read(path).unwrap(),
        b"{\"error_code\":\"structural_quality_abstention\",\"excluded_candidate_count\":0,\"fallback\":true,\"requested_selector\":\"factorized\",\"schema_version\":1,\"selector_used\":\"legacy\"}\n"
    );
}

#[test]
fn ineligible_preflight_rolls_back_without_calling_factorized_attempt() {
    let legacy = pipeline_choice_fixture();
    let decision = resolve_pipeline_choice(true, false, legacy.clone(), || {
        panic!("ineligible input must not construct features or call a model")
    });
    assert_eq!(decision.choice, legacy);
    assert!(matches!(decision.error, Some(FactorizedError::StructuralQualityAbstention)));
    assert_eq!(decision.excluded_candidate_count, 0);
}
```

Keep the existing success, flag-off, identity-fallback, and missing-Peeling tests.

- [ ] **Step 2: Run focused tests and capture RED**

```bash
cargo test -p sword2-lib sword::factorized_ranker::tests::structural_abstention_has_stable_error_code_and_status_bytes -- --nocapture
cargo test -p sword2-lib sword::tests::ineligible_preflight_rolls_back_without_calling_factorized_attempt -- --nocapture
```

Expected: failures because the error variant and eligibility argument are absent.

- [ ] **Step 3: Add the stable error variant**

Add to `FactorizedError`:

```rust
#[error("factorized selector abstained because required backbone atoms are incomplete")]
StructuralQualityAbstention,
```

and map only that variant to `"structural_quality_abstention"` in `code()`.

- [ ] **Step 4: Gate the pipeline transaction before the closure**

Change the pure decision seam to:

```rust
fn resolve_pipeline_choice<F>(
    use_factorized_ranker: bool,
    structurally_eligible: bool,
    legacy: PipelineChoice,
    factorized: F,
) -> PipelineDecision
where
    F: FnOnce() -> Result<(PipelineChoice, usize), factorized_ranker::FactorizedError>,
{
    if !use_factorized_ranker {
        return PipelineDecision { choice: legacy, error: None, excluded_candidate_count: 0 };
    }
    if !structurally_eligible {
        return PipelineDecision {
            choice: legacy,
            error: Some(factorized_ranker::FactorizedError::StructuralQualityAbstention),
            excluded_candidate_count: 0,
        };
    }
    // Existing transactional success/error match follows unchanged.
}
```

Parse the cleaned PDB once before DSSP, inspect its only processing chain, and carry the boolean through both the no-Peeling early return and normal decision path. If both structural ineligibility and missing Peeling are present, structural abstention takes precedence because model eligibility was already false. If structural evidence is complete, retain the current `feature_missing_context` result for missing Peeling/provenance.

Apply the same preflight before the development-only `SWORD2_DUMP_CANDIDATES` extraction block. An ineligible chain may emit its quality report and a stable dump rejection, but must not enter feature extraction even when the dump environment variable is set.

- [ ] **Step 5: Prove eligible behavior and error distinctions GREEN**

```bash
cargo test -p sword2-lib sword::tests::factorized_pipeline -- --nocapture
cargo test -p sword2-lib sword::factorized_ranker -- --nocapture
cargo test -p sword2-lib sword::tests:: -- --nocapture
```

Expected: no model closure call for ineligible fixtures; eligible success and later missing-context behavior remain unchanged.

- [ ] **Step 6: Commit Task 3**

```bash
git add sword2-lib/src/sword/factorized_ranker/mod.rs sword2-lib/src/sword/mod.rs
git diff --cached --check
git commit -m "fix: abstain on incomplete factorized structures"
```

## Task 4: Create the canonical input-only eligibility freeze

**Files:**
- Create: `benchmark/factorized_ranker/eligibility.py`
- Create: `benchmark/freeze_factorized_eligibility.py`
- Create: `benchmark/tests/test_factorized_eligibility.py`

- [ ] **Step 1: Write failing canonical-report and manifest tests**

The tests must define two synthetic identities and inject an inspector callback so no real SWORD or CATH input is used:

```python
def _report(*, eligible: bool) -> bytes:
    missing = [] if eligible else [{
        "author_residue_number": 7,
        "chain_id": "A",
        "missing_atoms": ["O"],
    }]
    return canonical_json_bytes({
        "candidate_residue_count": 10,
        "chain_id": "A",
        "complete_backbone_residue_count": 10 if eligible else 9,
        "eligible": eligible,
        "incomplete_residues": missing,
        "policy": "strict_complete_backbone_v1",
        "reason_code": None,
        "schema_version": 1,
        "structural_coverage": 1.0 if eligible else 0.9,
    })

def test_build_manifest_is_order_stable_and_partitions_all_ids(tmp_path):
    identities = [DatasetIdentity("1aaa", "oneA", "A"), DatasetIdentity("2bbb", "twoB", "B")]
    reports = {"oneA": _report(eligible=True), "twoB": _report(eligible=False)}
    manifest = build_eligibility_manifest(
        identities=identities,
        chain_root=tmp_path,
        dataset_sha256="d" * 64,
        runtime=synthetic_runtime(),
        runtime_manifest_sha256="c" * 64,
        binary_sha256="b" * 64,
        inspector=lambda identity, _path: reports[identity.entry_id],
        jobs=2,
    )
    assert manifest["dataset_ids"] == ["oneA", "twoB"]
    assert manifest["eligible_ids"] == ["oneA"]
    assert manifest["ineligible_ids"] == ["twoB"]
    assert manifest["eligible_count"] + manifest["ineligible_count"] == 2
    validate_eligibility_manifest(manifest)
```

Add one-field-at-a-time rejection tests for duplicate IDs, changed structure hash, changed quality bytes/hash, bad policy, bad runtime/binary hash, overlap/omission between eligible and ineligible sets, noncanonical JSON, a report chain mismatch, and inspector stderr/nonzero exit.

Add an identity-prefix test whose fourth field contains an unmistakable sentinel CATH label and assert the parser returns only `(pdb_id, entry_id, chain_id)`.

- [ ] **Step 2: Run tests and capture RED**

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_eligibility.py -q
```

Expected: collection/import error because the module and CLI do not exist.

- [ ] **Step 3: Implement strict canonical report validation**

Define exact report keys and checks:

```python
QUALITY_REPORT_KEYS = {
    "candidate_residue_count", "chain_id", "complete_backbone_residue_count",
    "eligible", "incomplete_residues", "policy", "reason_code",
    "schema_version", "structural_coverage",
}
POLICY = "strict_complete_backbone_v1"
REQUIRED_ATOM_ORDER = ("N", "CA", "C", "O")
```

`load_quality_record_bytes` must reject duplicate keys, nonfinite JSON constants, noncanonical serialization, unknown/missing keys, non-integer counts, coverage unequal to `complete/candidate` (or `0.0` for zero), malformed residue identities, missing-atom order differing from the `REQUIRED_ATOM_ORDER` subsequence, and any inconsistency among counts, `eligible`, incomplete rows, and `reason_code`.

- [ ] **Step 4: Implement identity-only metadata parsing**

Read stable bytes and parse only the first three comma-delimited byte fields:

```python
@dataclass(frozen=True)
class DatasetIdentity:
    pdb_id: str
    entry_id: str
    chain_id: str

def read_dataset_identity_prefixes(path: Path) -> tuple[list[DatasetIdentity], str]:
    data = stable_file_bytes(path)
    identities = []
    for line in data.splitlines():
        if not line or line.startswith(b"#"):
            continue
        prefix = line.split(b",", 3)
        if len(prefix) != 4:
            raise ValueError("dataset identity row is truncated")
        pdb_id, entry_id, chain_id = (field.decode("ascii") for field in prefix[:3])
        identities.append(validate_identity(DatasetIdentity(pdb_id, entry_id, chain_id)))
    return identities, hashlib.sha256(data).hexdigest()
```

Never pass the fourth remainder to `csv.reader`, `CathEntry`, `load_dataset`, or a log.

- [ ] **Step 5: Implement deterministic parallel inspection and manifest validation**

Production inspection runs:

```text
<binary> inspect-factorized-eligibility --input <cache>/chains/<entry_id>.pdb --chain <chain_id> --nmr-model 1
```

Use `ThreadPoolExecutor(max_workers=jobs)` only to wait for independent subprocesses; store results in source dataset order. Require return code zero, empty stderr, canonical stdout, report chain equal to metadata chain, unchanged input hash, and unchanged binary/runtime hashes before and after all workers.

The canonical manifest contains exactly:

```python
ELIGIBILITY_MANIFEST_KEYS = {
    "binary_sha256", "dataset", "dataset_id_count", "dataset_id_set_sha256",
    "dataset_ids", "dataset_sha256", "eligible_count", "eligible_id_set_sha256",
    "eligible_ids", "ineligibility_reason_counts", "ineligible_count",
    "ineligible_id_set_sha256", "ineligible_ids", "policy",
    "quality_record_sha256s", "quality_record_tree_sha256", "quality_records",
    "runtime_manifest_sha256", "runtime_source_git_commit", "schema_version",
    "structure_sha256s", "structure_tree_sha256",
}
```

Each `quality_record_sha256` hashes the standalone canonical report bytes including its final LF. The manifest itself is timestamp-free, canonical, written absent-only, and validated again after installation.

- [ ] **Step 6: Implement `create` and `verify` CLI commands**

The exact interface is:

```text
python -m benchmark.freeze_factorized_eligibility create \
  --dataset-metadata PATH --cache-dir PATH --runtime-manifest PATH \
  --binary PATH --repo-root PATH --jobs N --out PATH

python -m benchmark.freeze_factorized_eligibility verify \
  --manifest PATH --dataset-metadata PATH --cache-dir PATH \
  --runtime-manifest PATH --binary PATH --repo-root PATH
```

Both print only the lowercase manifest SHA-256 on stdout. Invalid evidence prints one concise stderr line and exits 2. Creation refuses an existing/symlink output and any output/input alias.

- [ ] **Step 7: Run focused tests GREEN**

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_eligibility.py -q
benchmark/.venv/bin/python -m py_compile \
  benchmark/factorized_ranker/eligibility.py \
  benchmark/freeze_factorized_eligibility.py
```

Expected: all focused tests pass.

- [ ] **Step 8: Commit Task 4**

```bash
git add benchmark/factorized_ranker/eligibility.py \
  benchmark/freeze_factorized_eligibility.py \
  benchmark/tests/test_factorized_eligibility.py
git diff --cached --check
git commit -m "feat: freeze factorized structural eligibility"
```

## Task 5: Evolve the runtime freeze without weakening old evidence parsing

**Files:**
- Modify: `benchmark/factorized_ranker/runtime_freeze.py:20-115,740-1035`
- Modify: `benchmark/tests/test_factorized_runtime_freeze.py`

- [ ] **Step 1: Write failing runtime-v2 tests**

Add tests asserting:

```python
assert created["schema_version"] == 2
assert set(created["evidence_tool_sha256s"]) == set(EVIDENCE_TOOL_PATHS_V2)
assert "benchmark/factorized_ranker/eligibility.py" in EVIDENCE_TOOL_PATHS_V2
assert "benchmark/freeze_factorized_eligibility.py" in EVIDENCE_TOOL_PATHS_V2
```

Retain a synthetic schema-v1 fixture and assert its schema is parsed against the v1 path set rather than silently reinterpreted as v2. Add tamper tests for either new tool hash and for changing the schema number without changing the path set.

- [ ] **Step 2: Run the focused test and capture RED**

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_runtime_freeze.py -q
```

Expected: failures because runtime schema v2/path sets are absent.

- [ ] **Step 3: Implement versioned evidence path sets**

Rename the current tuple to `EVIDENCE_TOOL_PATHS_V1`, define:

```python
EVIDENCE_TOOL_PATHS_V2 = tuple(sorted({
    *EVIDENCE_TOOL_PATHS_V1,
    "benchmark/factorized_ranker/eligibility.py",
    "benchmark/freeze_factorized_eligibility.py",
}))
RUNTIME_FREEZE_SCHEMA_VERSION = 2
```

Creation always emits schema 2. Validation selects the exact path set from the schema and rejects any other integer. Do not add model artifacts or weights to the new path set; they remain bound through the unchanged model manifest.

- [ ] **Step 4: Run runtime and eligibility suites GREEN**

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_runtime_freeze.py \
  benchmark/tests/test_factorized_eligibility.py -q
```

- [ ] **Step 5: Commit Task 5**

```bash
git add benchmark/factorized_ranker/runtime_freeze.py \
  benchmark/tests/test_factorized_runtime_freeze.py
git diff --cached --check
git commit -m "feat: bind eligibility tooling in runtime freeze"
```

## Task 6: Make locked collection accept only predeclared abstentions

**Files:**
- Modify: `benchmark/run_benchmark.py:55-170,245-265,889-1170,1340-1455,1768-2080`
- Modify: `benchmark/tests/test_run_benchmark.py`

- [ ] **Step 1: Write failing parser and run-policy tests**

Add tests for these exact cases:

```python
def test_structural_abstention_status_is_valid_only_when_success_not_required():
    data = canonical_json_bytes({
        "error_code": "structural_quality_abstention",
        "excluded_candidate_count": 0,
        "fallback": True,
        "requested_selector": "factorized",
        "schema_version": 1,
        "selector_used": "legacy",
    })
    parsed = parse_selector_status(data, "factorized", require_success=False)
    assert parsed["error_code"] == "structural_quality_abstention"
    with pytest.raises(ValueError):
        parse_selector_status(data, "factorized", require_success=True)

def _factorized_row(entry_id: str, *, success: bool, code: str = "structural_quality_abstention", exclusions: int = 0):
    return SimpleNamespace(
        entry_id=entry_id,
        requested_selector="factorized",
        selector_used="factorized" if success else "legacy",
        fallback=not success,
        error_code="" if success else code,
        selector_warning_code="" if success else "factorized_fallback",
        excluded_candidate_count=exclusions,
    )

def test_factorized_outcome_scores_eligible_and_retains_ineligible_run_only():
    eligible = frozenset({"a"})
    assert _validate_locked_factorized_outcome("a", _factorized_row("a", success=True), eligible)
    assert not _validate_locked_factorized_outcome("b", _factorized_row("b", success=False), eligible)

@pytest.mark.parametrize(
    ("entry_id", "row"),
    [
        ("a", _factorized_row("a", success=False)),
        ("b", _factorized_row("b", success=True)),
        ("b", _factorized_row("b", success=False, code="feature_missing_context")),
        ("b", _factorized_row("b", success=False, exclusions=1)),
    ],
)
def test_expected_abstention_rejects_wrong_id_error_or_exclusion(entry_id, row):
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_factorized_outcome(entry_id, row, frozenset({"a"}))

def test_resource_assignments_use_only_eligible_ids():
    assignments = _eligible_resource_assignments(
        dataset_ids=("a", "b", "c"),
        eligible_ids=frozenset({"a", "c"}),
    )
    assert set(assignments) == {"a", "c"}
    assert "b" not in assignments

def test_locked_worker_contract_is_role_specific():
    assert _validate_locked_jobs("factorized-accuracy", 32) == 32
    assert _validate_locked_jobs("legacy-accuracy", 32) == 32
    assert _validate_locked_jobs("paired-sword-resources", 1) == 1
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_jobs("factorized-accuracy", 31)
    with pytest.raises(LockedBenchmarkError):
        _validate_locked_jobs("paired-sword-resources", 2)
```

Also require `--locked-eligibility-manifest` for all three locked roles and `--locked-jobs 32` for accuracy versus `--locked-jobs 1` for resources.

- [ ] **Step 2: Run focused tests and capture RED**

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_run_benchmark.py -q
```

Expected: failures for missing arguments, preflight fields, and abstention policy.

- [ ] **Step 3: Extend locked preflight authority**

Add CLI arguments:

```python
parser.add_argument("--locked-eligibility-manifest", type=Path, default=None)
parser.add_argument("--locked-jobs", type=int, default=None)
```

Extend `LockedPreflight` with the canonical eligibility payload/hash and exact eligible/ineligible tuples. Call `verify_eligibility_manifest` using the already verified runtime, binary, dataset, cache, and repository paths. Require its complete identity, structure-hash, runtime-hash, and binary-hash graph to equal preflight's independently derived graph.

Require `locked_jobs == 32` for `legacy-accuracy` and `factorized-accuracy`; require `locked_jobs == 1` for `paired-sword-resources`. Preserve `sword2_threads == 1` and all existing BLAS/Rayon environment constraints.

- [ ] **Step 4: Add exact status helpers**

Keep `_require_locked_sword_success` unchanged and add:

```python
def _require_locked_structural_abstention(row: LockedRunRow) -> None:
    if not (
        row.requested_selector == "factorized"
        and row.selector_used == "legacy"
        and row.fallback is True
        and row.error_code == "structural_quality_abstention"
        and row.selector_warning_code == "factorized_fallback"
        and row.excluded_candidate_count == 0
    ):
        raise LockedBenchmarkError("locked structural abstention is inconsistent")
```

Do not infer abstention from logs or `residue_quality.json`; selector status and the frozen eligibility ID set must agree.

- [ ] **Step 5: Update factorized accuracy collection**

Run all 663 factorized requests with a 32-worker `ThreadPoolExecutor`. Each worker writes only its private raw directory and returns predictions plus one `LockedRunRow`; the main thread performs status validation, scoring, and atomic evidence CSV updates.

For an eligible ID: require actual factorized success and exactly one rank-one score path. For an ineligible ID: require exact structural abstention, retain its run and raw legacy summary for audit, and emit no `ScoreRow`. Any other combination raises and preserves an incomplete manifest-free result directory.

On the legacy role, require successful legacy SWORD for all 663 before competitor collection; its metric rows remain full-population evidence.

- [ ] **Step 6: Update resource collection and manifest schema**

Build counterbalanced assignments from `eligible_ids` only and keep the existing nested sequential loop. Both selectors must succeed for every resource-pair ID. Add to locked manifest schema v2:

```python
"eligibility_manifest_sha256",
"eligibility_policy",
"factorized_eligible_count",
"factorized_eligible_id_set_sha256",
"structural_abstention_count",
"structural_abstention_id_set_sha256",
"locked_jobs",
```

For factorized accuracy, `fallback_count == structural_abstention_count`; for legacy and resources it is zero. `selector_counts` continues to count requested commands, while validators inspect `selector_used` on every row. Resource success counts equal the eligible set, not 663.

- [ ] **Step 7: Run runner tests GREEN**

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_run_benchmark.py \
  benchmark/tests/test_runners.py -q
```

Expected: expected abstentions continue, unexpected fallbacks stop, resource ordering is unchanged within the eligible set, and concurrency is bounded.

- [ ] **Step 8: Commit Task 6**

```bash
git add benchmark/run_benchmark.py benchmark/tests/test_run_benchmark.py
git diff --cached --check
git commit -m "feat: collect predeclared factorized abstentions"
```

## Task 7: Validate coverage before conditional metrics

**Files:**
- Modify: `benchmark/evaluate_factorized_acceptance.py:35-180,330-535,600-1435,1840-2165`
- Modify: `benchmark/tests/test_factorized_acceptance.py`

- [ ] **Step 1: Convert the synthetic evidence fixture to eligible/ineligible sets**

Use a three-ID fixture with `eligible={a,c}`, `ineligible={b}`. Factorized runs contain success rows for `a,c` and exact abstention for `b`; factorized scores contain only `a,c`; legacy/Merizo contain all three; Chainsaw expected success contains `a,b`; resources contain both selectors for `a,c` only.

Assert coverage output contains no metric names and exactly:

```python
assert payload["dataset_id_count"] == 3
assert payload["factorized_eligible_count"] == 2
assert payload["structural_abstention_count"] == 1
assert payload["factorized_structural_coverage"] == 2 / 3
assert payload["resource_pair_count"] == 2
assert payload["fallback_count"] == 1
```

- [ ] **Step 2: Add failing tamper tests**

Parameterize rejection of:

- score row for ineligible `b`;
- missing score for eligible `a`;
- success status for ineligible `b`;
- abstention status for eligible `a`;
- wrong abstention error/exclusion count;
- resource row for `b` or missing resource pair for `c`;
- eligibility manifest hash/policy/input hash mutation;
- overlap, omission, or reordering in eligibility IDs;
- conditional Merizo/standalone ID mismatch;
- metric access before coverage succeeds.

- [ ] **Step 3: Run focused tests and capture RED**

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests/test_factorized_acceptance.py -q
```

Expected: current schema requires full factorized/resource score coverage and zero fallback.

- [ ] **Step 4: Implement locked manifest and coverage schema v2**

Add `eligibility_manifest` to `EVIDENCE_ARGUMENT_NAMES` and validate it before score metrics. Require all three role manifests to bind its exact hash/policy/sets. Role coverage becomes:

```text
legacy accuracy:       successful SWORD + Merizo on all 663; Chainsaw on frozen success set
factorized accuracy:   663 run rows; actual factorized scores on eligible set; exact abstentions on complement
paired resources:      one successful legacy/factorized pair for every eligible ID only
```

Coverage schema v2 adds:

```python
"eligibility_manifest_sha256",
"eligibility_policy",
"factorized_eligible_count",
"factorized_eligible_id_set_sha256",
"structural_abstention_count",
"structural_abstention_id_set_sha256",
"factorized_structural_coverage",
```

Require `eligible + abstention == dataset_id_count`, `resource_pair_count == eligible`, and `fallback_count == abstention`. Do not require fallback count zero.

- [ ] **Step 5: Restrict every conditional comparison to the same eligible set**

After revalidating coverage, select:

```python
eligible_ids = set(eligibility["eligible_ids"])
factorized = select_factorized_scores(expected_ids=eligible_ids)
legacy = select_legacy_scores(expected_ids=all_ids).loc[sorted(eligible_ids)]
merizo = select_merizo_scores(expected_ids=all_ids).loc[sorted(eligible_ids)]
chainsaw_ids = frozen_chainsaw_ids & eligible_ids
chainsaw = select_chainsaw_scores(expected_ids=frozen_chainsaw_ids).loc[sorted(chainsaw_ids)]
standalone = load_standalone_metrics(all_ids).loc[sorted(eligible_ids)]
resources = resource_metric_frames(expected_ids=eligible_ids)
```

Truth continuity cohorts are loaded only after coverage and then restricted to eligible IDs. The nine existing thresholds, comparisons, bootstrap method, seed, and replicate count stay unchanged, but every gate denominator is labeled conditional.

- [ ] **Step 6: Add transparent coverage and population diagnostics**

Acceptance schema v2 retains the nine gates and adds diagnostics for:

- full denominator, eligible count, abstention count, and structural coverage;
- ineligibility reason counts from the input-only manifest;
- conditional mean NDO for factorized, runtime legacy, Merizo, and available Chainsaw;
- full-population mean NDO for runtime legacy and Merizo, plus Chainsaw's frozen-success denominator;
- explicit `default_promotion_compatible = false` when structural coverage is below 1.0 or the runtime cache contract is incompatible.

Render structural coverage and conditional-denominator caveats above the gate table in Markdown. Never render an abstained chain's legacy fallback as a factorized score.

- [ ] **Step 7: Keep promotion fail-closed**

`validate-promotion` must rederive eligibility, coverage, conditional metrics, Markdown, all hashes, and the committed evidence graph. Even if all nine conditional gates pass, return `KEEP_OPT_IN` when coverage is below 100% or the runtime cache contract is non-promotable. Invalid evidence remains exit 2 and produces no decision token.

- [ ] **Step 8: Run evaluator and integration tests GREEN**

```bash
benchmark/.venv/bin/python -m pytest \
  benchmark/tests/test_factorized_acceptance.py \
  benchmark/tests/test_run_benchmark.py \
  benchmark/tests/test_factorized_runtime_freeze.py \
  benchmark/tests/test_factorized_eligibility.py \
  benchmark/tests/test_runners.py -q
```

- [ ] **Step 9: Commit Task 7**

```bash
git add benchmark/evaluate_factorized_acceptance.py \
  benchmark/tests/test_factorized_acceptance.py
git diff --cached --check
git commit -m "feat: evaluate eligible factorized coverage"
```

## Task 8: Complete development-only verification

**Files:**
- Verify only: all source/test files above
- Create locally, do not stage: `.superpowers/sdd/2026-08-07-factorized-structural-ranker/task-19-report.md`

- [ ] **Step 1: Run scoped formatting without touching unrelated files**

```bash
cargo fmt --check --package sword2-lib --package sword
```

If the repository-wide formatter exposes pre-existing unrelated drift, format only changed Rust files with the pinned rustfmt command, then rerun `git diff --check`. Never stage unrelated formatting.

- [ ] **Step 2: Run Rust verification**

Run independent compile/test gates with available cores where Cargo safely supports it:

```bash
cargo check --workspace
cargo test --workspace
cargo build --locked --release --bin sword2
```

Expected: CLI tests, library tests, DSSP golden, and factorized golden all pass; only documented pre-existing warnings may remain.

- [ ] **Step 3: Run full Python verification**

```bash
benchmark/.venv/bin/python -m pytest benchmark/tests -q
```

Expected: full suite passes.

- [ ] **Step 4: Run development-only behavior smokes**

Use fresh temporary directories and CATH-17287 cached structures only:

- complete example `12e8H`: inspector eligible; factorized status succeeds;
- missing-O example `1a5tA`: inspector ineligible; factorized status is exact structural abstention; legacy summary still exists;
- missing-CA example `1bgwA`: same exact abstention;
- empty-Peeling example `2wjvD`: if structurally complete, retain `feature_missing_context`; if structurally incomplete, structural abstention takes precedence as specified.

For each normal run, compare `residue_quality.json` bytes with the read-only inspector bytes. Confirm no learned feature dump/model status is produced for abstained inputs.

- [ ] **Step 5: Recompute the CATH-17287 input-only profile**

Run the new eligibility creator with `--jobs 32` into a temporary absent path. Require exactly 17,286 identities for the current hash-frozen external metadata (the historical `CATH-17287` name is not its row count), exact accounting, deterministic byte identity across a second location, and reconcile against the prior 591-chain incomplete diagnostic without reading any CATH-663 input. Record counts/hashes in the local report and explain any diagnostic-count difference rather than forcing the stale nominal count.

- [ ] **Step 6: Remove only temporary development artifacts**

Resolve and validate each exact temporary path under `/tmp` before deleting it with `find <exact-path> -depth -delete`. Preserve committed source, model artifacts, old failed ledgers, and unrelated user files.

- [ ] **Step 7: Record final source state**

```bash
git diff --check
git status --short
git log --oneline -10
```

If verification required a source/test correction, use RED/GREEN and a surgical fix commit before continuing. The source commit used by runtime freeze must have clean tracked runtime/evidence paths.

## Task 9: Freeze runtime and eligibility before replacement evaluation

**Files:**
- Create and commit: `benchmark/models/factorized_ranker_v1_runtime_v2.json`
- Create and commit before CATH inspection: `benchmark/models/factorized_ranker_v1_eligibility_intent_attempt3.json`
- Create and commit after input-only inspection: `benchmark/models/factorized_ranker_v1_structural_eligibility.json`
- Create and commit before model/competitor execution: `benchmark/models/factorized_ranker_v1_locked_intent_attempt3.json`
- Preserve read-only: all attempt-1/attempt-2 ledgers and intents

- [ ] **Step 1: Create and verify the runtime-v2 manifest**

With clean runtime/evidence paths and the just-built release binary:

```bash
benchmark/.venv/bin/python -m benchmark.freeze_factorized_runtime create \
  --model-manifest benchmark/models/factorized_ranker_v1_manifest.json \
  --corpus-dir benchmark/data/cath17287_factorized_corpus_v1 \
  --fold-manifest benchmark/models/cath17287_factorized_folds_v1.json \
  --oof-predictions benchmark/models/factorized_ranker_v1_oof.csv \
  --binary target/release/sword2 \
  --repo-root . \
  --out benchmark/models/factorized_ranker_v1_runtime_v2.json

benchmark/.venv/bin/python -m benchmark.freeze_factorized_runtime verify \
  --manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --binary target/release/sword2 --repo-root .
```

Commit only the new runtime manifest with message `model: refreeze factorized runtime with structural abstention`.

- [ ] **Step 2: Commit a one-shot eligibility-inspection intent**

Create canonical `factorized_ranker_v1_eligibility_intent_attempt3.json` containing schema/version, `authorized_once`, dataset name and already-known metadata/ID/input-tree hashes, runtime/model/binary/source hashes, policy, exact normalized inspector command, `jobs=32`, absent output path, the prohibition on labels/scores/competitor outputs, and the prior-attempt audit disclosure. It contains no eligibility count.

Validate canonical bytes and commit only this intent with message `bench: authorize structural eligibility inspection`.

- [ ] **Step 3: Perform the only new pre-evaluation CATH-663 operation**

Run:

```bash
benchmark/.venv/bin/python -m benchmark.freeze_factorized_eligibility create \
  --dataset-metadata /home/chili/cretin/PROJECTS/Merizo/datasets/merizo_domains/CATH-663.csv \
  --cache-dir benchmark/cache \
  --runtime-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --binary target/release/sword2 \
  --repo-root . --jobs 32 \
  --out benchmark/models/factorized_ranker_v1_structural_eligibility.json
```

Then run `verify` with the same authority inputs. Do not open CATH labels, existing CATH score CSVs, factorized predictions, or competitor outputs. Report only eligible/ineligible counts and input-quality reasons.

Commit only the eligibility manifest with message `bench: freeze factorized structural eligibility`.

- [ ] **Step 4: Create the replacement locked intent**

Create canonical `factorized_ranker_v1_locked_intent_attempt3.json` bound to the runtime-v2 and eligibility manifests, old model/standalone/Chainsaw artifacts, exact external Merizo/Chainsaw hashes, full 663 denominator, eligible resource/gate denominator, expected abstention complement, 32-worker accuracy commands, one-worker paired-resource command, unchanged nine gates, no-timeout invalidation policy, new attempt-3 result paths, and disclosure that this is a protocol-amended replacement after status-level failure.

Commit only this intent with message `bench: authorize eligibility-aware factorized evaluation`.

- [ ] **Step 5: Reverify every frozen authority before launch**

Require runtime, eligibility, model, intent, dataset, structure tree, standalone baseline, expected Chainsaw set, external artifacts, binary, source, Cargo.lock, Python environments, CUDA identity, free disk/inodes, and 64-GB memory readiness to match the intent. Any mismatch stops without running a model.

## Task 10: Run the replacement locked protocol and retain opt-in

**Files:**
- Create: attempt-3 resource, legacy, factorized result directories
- Create: attempt-3 local append-only ledger and logs
- Create only after coverage succeeds: attempt-3 coverage and acceptance artifacts
- Update only after final validator succeeds: `benchmark/REPORT.md` and user-facing status documentation

- [ ] **Step 1: Start an inspectable append-only ledger**

Before commands, create a fresh attempt-3 ledger containing exact intent bytes/hash, environment snapshot, external artifact descriptors, and append-only `events.jsonl`. Each event records UTC time, role, phase, completed/total chains, process count, command hash, and status without metric values. Keep command stdout/stderr in named ledger files so progress can be checked while running.

- [ ] **Step 2: Run sequential counterbalanced resource pairs**

Set the frozen process environment once:

```bash
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RAYON_NUM_THREADS=1
```

Then run the exact resource command into an absent path:

```bash
benchmark/.venv/bin/python -m benchmark.run_benchmark \
  --dataset cath663 \
  --cache-dir benchmark/cache \
  --results-dir benchmark/results_factorized_locked_resources_attempt3 \
  --tools sword2-rust \
  --sword2-threads 1 \
  --locked-jobs 1 \
  --no-download --strict \
  --locked-factorized-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --locked-eligibility-manifest benchmark/models/factorized_ranker_v1_structural_eligibility.json \
  --locked-role paired-sword-resources \
  >benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/paired.stdout \
  2>benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/paired.stderr
```

It processes eligible IDs only, in the precommitted seed-37 counterbalanced order. Do not run CATH as a warm-up, flush OS caches, impose a timeout, or parallelize resource pairs.

- [ ] **Step 3: Run accuracy collection with bounded concurrency**

Run the legacy/competitor role and factorized role using 32 independent one-thread SWORD children:

```bash
benchmark/.venv/bin/python -m benchmark.run_benchmark \
  --dataset cath663 \
  --cache-dir benchmark/cache \
  --results-dir benchmark/results_factorized_locked_legacy_attempt3 \
  --tools sword2-rust,merizo,chainsaw \
  --sword2-threads 1 \
  --locked-jobs 32 \
  --no-download --strict \
  --locked-factorized-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --locked-eligibility-manifest benchmark/models/factorized_ranker_v1_structural_eligibility.json \
  --locked-role legacy-accuracy \
  --chainsaw-expected-success-ids benchmark/models/cath663_chainsaw_expected_success_ids_v1.txt \
  --locked-artifact merizo_python=/home/chili/cretin/PROJECTS/Merizo/merizo/bin/python \
  --locked-artifact merizo_source_tree=/home/chili/cretin/PROJECTS/Merizo \
  --locked-artifact merizo_model_tree=/home/chili/cretin/PROJECTS/Merizo/model \
  --locked-artifact chainsaw_python=/home/chili/cretin/PROJECTS/chainsaw/chswEnv/bin/python \
  --locked-artifact chainsaw_source_tree=/home/chili/cretin/PROJECTS/chainsaw \
  --locked-artifact chainsaw_model_tree=/home/chili/cretin/PROJECTS/chainsaw/saved_models/model_v3 \
  >benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/legacy.stdout \
  2>benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/legacy.stderr

benchmark/.venv/bin/python -m benchmark.run_benchmark \
  --dataset cath663 \
  --cache-dir benchmark/cache \
  --results-dir benchmark/results_factorized_locked_model_attempt3 \
  --tools sword2-rust \
  --sword2-threads 1 \
  --locked-jobs 32 \
  --sword2-extra-args=--use-factorized-ranker \
  --no-download --strict \
  --locked-factorized-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --locked-eligibility-manifest benchmark/models/factorized_ranker_v1_structural_eligibility.json \
  --locked-role factorized-accuracy \
  >benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/factorized.stdout \
  2>benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/factorized.stderr
```

Merizo and Chainsaw use the previously hashed installations at `../Merizo` and `../chainsaw`. Expected structural abstentions continue; any unexpected fallback or set mismatch invalidates the attempt.

- [ ] **Step 4: Validate coverage without reading metrics**

Run the exact coverage-only command:

```bash
benchmark/.venv/bin/python -m benchmark.evaluate_factorized_acceptance coverage \
  --legacy-manifest benchmark/results_factorized_locked_legacy_attempt3/benchmark_manifest.json \
  --legacy-scores benchmark/results_factorized_locked_legacy_attempt3/scores.csv \
  --legacy-runs benchmark/results_factorized_locked_legacy_attempt3/runs.csv \
  --legacy-failures benchmark/results_factorized_locked_legacy_attempt3/failures.csv \
  --factorized-manifest benchmark/results_factorized_locked_model_attempt3/benchmark_manifest.json \
  --factorized-scores benchmark/results_factorized_locked_model_attempt3/scores.csv \
  --factorized-runs benchmark/results_factorized_locked_model_attempt3/runs.csv \
  --factorized-failures benchmark/results_factorized_locked_model_attempt3/failures.csv \
  --resource-manifest benchmark/results_factorized_locked_resources_attempt3/benchmark_manifest.json \
  --resource-runs benchmark/results_factorized_locked_resources_attempt3/runs.csv \
  --resource-failures benchmark/results_factorized_locked_resources_attempt3/failures.csv \
  --dataset-metadata /home/chili/cretin/PROJECTS/Merizo/datasets/merizo_domains/CATH-663.csv \
  --chainsaw-expected-success-ids benchmark/models/cath663_chainsaw_expected_success_ids_v1.txt \
  --model-manifest benchmark/models/factorized_ranker_v1_manifest.json \
  --runtime-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --eligibility-manifest benchmark/models/factorized_ranker_v1_structural_eligibility.json \
  --standalone-baseline benchmark/models/cath663_standalone_structural_baseline.csv \
  --coverage-out benchmark/models/factorized_ranker_v1_coverage_attempt3.json \
  >benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/coverage.stdout \
  2>benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/coverage.stderr
```

It must prove full 663 run identity, exact eligible/ineligible statuses, eligible-only factorized score identities, eligible-only resource pairs, full legacy/Merizo identities, frozen Chainsaw identities, hashes, environments, commands, and manifests before writing a canonical coverage attestation.

If coverage exits 2, preserve all bytes, append invalidation evidence, stop, and do not evaluate metrics.

- [ ] **Step 5: Evaluate the nine conditional gates once**

Only after coverage exit 0, run:

```bash
benchmark/.venv/bin/python -m benchmark.evaluate_factorized_acceptance evaluate \
  --coverage-attestation benchmark/models/factorized_ranker_v1_coverage_attempt3.json \
  --legacy-manifest benchmark/results_factorized_locked_legacy_attempt3/benchmark_manifest.json \
  --legacy-scores benchmark/results_factorized_locked_legacy_attempt3/scores.csv \
  --legacy-runs benchmark/results_factorized_locked_legacy_attempt3/runs.csv \
  --legacy-failures benchmark/results_factorized_locked_legacy_attempt3/failures.csv \
  --factorized-manifest benchmark/results_factorized_locked_model_attempt3/benchmark_manifest.json \
  --factorized-scores benchmark/results_factorized_locked_model_attempt3/scores.csv \
  --factorized-runs benchmark/results_factorized_locked_model_attempt3/runs.csv \
  --factorized-failures benchmark/results_factorized_locked_model_attempt3/failures.csv \
  --resource-manifest benchmark/results_factorized_locked_resources_attempt3/benchmark_manifest.json \
  --resource-runs benchmark/results_factorized_locked_resources_attempt3/runs.csv \
  --resource-failures benchmark/results_factorized_locked_resources_attempt3/failures.csv \
  --dataset-metadata /home/chili/cretin/PROJECTS/Merizo/datasets/merizo_domains/CATH-663.csv \
  --chainsaw-expected-success-ids benchmark/models/cath663_chainsaw_expected_success_ids_v1.txt \
  --model-manifest benchmark/models/factorized_ranker_v1_manifest.json \
  --runtime-manifest benchmark/models/factorized_ranker_v1_runtime_v2.json \
  --eligibility-manifest benchmark/models/factorized_ranker_v1_structural_eligibility.json \
  --standalone-baseline benchmark/models/cath663_standalone_structural_baseline.csv \
  --json-out benchmark/models/factorized_ranker_v1_acceptance_attempt3.json \
  --markdown-out benchmark/FACTORIZED_RANKER_ACCEPTANCE_ATTEMPT3.md \
  --bootstrap-replicates 10000 --seed 37 \
  >benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/evaluate.stdout \
  2>benchmark/data/factorized_ranker_v1_locked_run_ledger_attempt3/evaluate.stderr
```

Store canonical JSON and Markdown. Report structural coverage separately from every conditional metric and include the full-population legacy/Merizo diagnostics.

- [ ] **Step 6: Run the complete promotion validator**

Revalidate committed intent/runtime/eligibility/model/input/output hashes, coverage, acceptance, Markdown, gate formulas, denominators, and source commit. The expected rollout decision is `KEEP_OPT_IN` whenever structural coverage is below 100% or the runtime cache contract remains non-promotable, even if all nine conditional gates pass.

- [ ] **Step 7: Commit final valid evidence and report**

Only valid, complete evidence is staged. Preserve prior failed ledgers. Append a report section stating:

- exact full denominator and structural coverage;
- conditional eligible-set results and denominators;
- abstention reasons;
- prior protocol amendment disclosure;
- nine gate outcomes;
- resource ratios;
- final `KEEP_OPT_IN` decision and prerequisites for a future gap-aware model.

Run final `git diff --check`, focused/full Rust and Python verification, then commit with message `bench: record eligibility-aware factorized evaluation`.

## Completion checklist

- [ ] Existing frozen weights and feature schema are byte-identical.
- [ ] Ordinary legacy selection behavior is unchanged.
- [ ] Inspector and normal-run quality bytes are identical.
- [ ] Ineligible input cannot call factorized feature extraction or either model head.
- [ ] Eligibility is frozen from input bytes before labels, predictions, or metrics.
- [ ] All 663 identities remain in coverage accounting.
- [ ] Factorized/competitor comparisons use the identical frozen eligible subset.
- [ ] Expected abstentions are never scored as factorized predictions.
- [ ] Unexpected fallbacks remain fail-closed.
- [ ] Accuracy collection uses 32 one-thread workers; resource collection is sequential.
- [ ] Full verification passes before every completion claim.
- [ ] Temporary development artifacts are deleted; audit ledgers are preserved.
- [ ] Final rollout remains opt-in unless a future model/protocol independently resolves both missing-structure coverage and cache determinism.
