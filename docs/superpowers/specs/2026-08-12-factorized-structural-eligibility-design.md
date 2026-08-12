# Factorized Structural Eligibility and Scientific Abstention

Date: 2026-08-12
Status: approved conversational design; written specification awaiting user review

## Context

The frozen factorized structural ranker requires one residue population shared by
CA geometry, DSSP, Peeling, candidate identity, and output numbering. The current
pipeline does not provide that invariant:

- chain cleaning retains every standard, non-insertion-code ATOM residue;
- CA geometry retains only residues containing CA;
- DSSP retains only residues containing N, CA, C, and O.

The first locked CATH-663 factorized attempt encountered this mismatch and
correctly stopped before metrics were evaluated. The structured status reported
`feature_missing_context`, then returned the legacy selection as a user-facing
fallback. That output is not factorized scientific evidence.

A development-only audit of the 591 rejected CATH-17287 chain dumps found 753
incomplete-backbone residues. Of those chains, 511 (86.5%) had one incomplete
residue, 561 (94.9%) had at most two, and 563 (95.3%) had at most 1% incomplete
residues. There were also material outliers: up to 20 residues and 6.1% of a
chain. Sparse missingness is therefore common enough that silently deleting all
affected chains would bias coverage, while silently filtering their residues
would apply the frozen model to an unvalidated representation.

## Decision

Version 1 of the frozen factorized ranker will use strict structural eligibility.
It will not delete incomplete residues and will not rank a chain containing an
eligible SWORD residue without complete N, CA, C, and O evidence. Such a chain
remains a valid SWORD2 input, may receive the legacy result for usability, and is
recorded as an explicit factorized abstention.

This is a model-validity decision, not a claim that the entire biological
structure is useless. Gap-tolerant learned inference requires a separately
trained and calibrated model version.

## Goals

1. Give the frozen model only the structural representation on which its feature
   contract was trained and validated.
2. Distinguish an actual factorized selection from a legacy fallback or
   structural-quality abstention in machine-readable output.
3. Keep every CATH-663 identity in the audit denominator and prevent post-hoc
   exclusion from improving reported metrics.
4. Determine eligibility solely from input structure bytes before reading labels,
   predictions, model scores, or benchmark metrics.
5. Preserve the ordinary legacy pipeline and current frozen model weights.
6. Produce enough residue-level evidence to reproduce every eligibility decision.

## Non-goals

- Do not impute missing coordinates or DSSP values.
- Do not remove incomplete residues and feed the compacted chain to the current
  model.
- Do not add missingness features, retrain trees, or change the 154-feature schema.
- Do not tune an eligibility threshold on CATH-663.
- Do not call a deterministic or legacy output a learned factorized decision.
- Do not delete or reuse the failed locked-attempt ledger as valid evidence.
- Do not make the factorized selector default-on in this change.

## Alternatives considered

### 1. Strict abstention for the frozen model

This design is selected. It preserves the meaning of the existing features and
weights, has a deterministic fail-closed contract, and avoids fabricating model
coverage. Its cost is less than 100% learned-selector coverage.

### 2. Complete-backbone filtering before current-model inference

This would make the internal arrays align and recover many chains. It is rejected
for version 1 because removing a residue changes contact probabilities, DSSP
context, Peeling candidates, domain-size features, and numbering. The frozen
trees were not calibrated for this intervention, and endpoint-based output
mapping can implicitly assign an omitted residue despite having no evidence for
it.

### 3. Gap-aware model and explicit missing-residue mask

This is the preferred future coverage extension. It would use development-only
synthetic deletion experiments, exact observed-residue mappings, missingness
features, and new frozen weights. It is not part of this repair because it is a
new model rather than a correction to the current model's execution contract.

## Authoritative structural eligibility contract

### Candidate residues

The inspector evaluates the same residues that `clean_chain_for_sword` currently
admits:

- the residue has at least one non-HETATM ATOM record;
- its residue name is one of the 20 standard amino acids;
- its insertion code is blank.

The change must not broaden or narrow that base residue set. Nonstandard,
HETATM-only, and insertion-code residues remain outside both the numerator and
denominator of the structural coverage calculation.

### Required evidence

Every candidate residue must have finite coordinates for N, CA, C, and O under
the same alternate-location acceptance used by the DSSP extractor: blank or
`A`. An atom present only in another alternate location does not satisfy version
1 eligibility. At least one accepted atom of each required name is sufficient;
unrelated side-chain atoms do not affect eligibility.

The implementation must expose one Rust predicate for this contract and reuse it
from inspection and factorized preflight. Tests must keep it in exact agreement
with DSSP's accepted atom semantics. A later alternate-conformer redesign would
require a new policy version.

### Chain decision

A chain is factorized-eligible exactly when:

- it has at least one candidate residue;
- every candidate residue has all four required atoms with finite coordinates.

Version 1 deliberately has no percentage threshold. One incomplete residue is
enough to abstain. This avoids selecting a cutoff without calibration.

An eligible chain can still abstain later for a different, explicitly named
reason such as missing Peeling iterations, missing provenance, or an unrankable
candidate lattice. Structural eligibility does not weaken the existing context
or finite-value checks.

## Quality record

The Rust inspector produces a canonical, schema-versioned record with:

- `schema_version = 1`;
- `policy = "strict_complete_backbone_v1"`;
- chain identity;
- candidate residue count;
- complete-backbone residue count;
- structural coverage as complete count divided by candidate count;
- `eligible`;
- one entry per incomplete residue containing original author residue number,
  chain identifier, and a `missing_atoms` subsequence in the canonical backbone
  order `N`, `CA`, `C`, `O`;
- a stable reason code for an empty candidate-residue population.

Structural coverage is `0.0` when the candidate-residue count is zero; otherwise
it is the exact complete count divided by the candidate count.

Records use canonical JSON: UTF-8, sorted object keys, no duplicate keys, finite
numbers only, and one trailing newline. Residue entries retain source order.
The report contains no model output, CATH label, benchmark truth, or metric.

Normal runs write `residue_quality.json` inside the chain result directory. The
read-only command

```text
sword2 inspect-factorized-eligibility --input PATH [--chain ID] [--nmr-model N]
```

emits the identical canonical record to stdout and performs no filesystem write.
It exits zero for both eligible and ineligible successfully parsed chains; input,
model, or chain-selection errors exit nonzero and emit no JSON. The record's
chain identity is the selected structure chain ID; benchmark entry identity is
bound externally in the eligibility manifest. The command does not run DSSP,
Peeling, candidate generation, energy scoring, or selection. The locked
benchmark uses this inspection path so its eligibility set is generated by the
same frozen Rust binary as runtime enforcement.

## Runtime data flow

1. Parse the requested chain and apply the existing base cleaning rules.
2. Build the structural-quality record before writing pipeline intermediates.
3. For an ordinary legacy request, preserve existing selection behavior. The
   quality report is informational and does not reject the legacy run.
4. For a factorized request on an ineligible chain, do not construct factorized
   features and do not invoke either embedded head.
5. Preserve the current user-facing legacy output, but write structured selector
   status with:
   - `requested_selector = "factorized"`;
   - `selector_used = "legacy"`;
   - `fallback = true`;
   - `error_code = "structural_quality_abstention"`;
   - zero factorized candidate exclusions, because ranking never began.
6. For an eligible chain, continue through the existing strict context checks and
   factorized selection transaction.
7. Any later context failure remains a fallback with its own existing error code;
   it must not be relabeled as structural ineligibility.

No learned probability, margin, count score, or winner score is emitted for an
abstained chain. Logs may summarize the missing atom count but are not the source
of scientific status.

## Empty or unrankable lattices

A chain with complete backbone evidence but no usable Peeling iterations or no
rankable candidate lattice is not an actual learned selection. Version 1 records
the existing context-specific abstention/fallback rather than declaring
factorized success. A future protocol may introduce a separately named
deterministic single-candidate outcome, but this repair does not conflate it with
model coverage.

## Locked benchmark protocol amendment

### Blinded eligibility freeze

After implementation, tests, commit, release build, and runtime freeze are
complete, create a replacement evaluation intent before inspecting further
CATH-663 inputs. Run only the frozen binary's read-only quality inspector over
the exact 663 canonical structure files. Freeze:

- ordered complete ID set and its hash;
- ordered eligible and ineligible ID sets and their hashes;
- per-input structure hash;
- per-ID canonical quality-record hash;
- policy, binary, source, model, and runtime-manifest hashes.

This phase may inspect atoms but must not load CATH chopping labels, competitor
predictions, model predictions, or metric columns. Eligibility is immutable once
the input-only manifest is installed.

### Collection

All 663 identities remain represented in locked run evidence.

- Legacy SWORD2 and competitors retain the full 663-chain collection contract.
- Factorized collection attempts every ID.
- Every frozen eligible ID must produce an actual factorized status and score row.
- Every frozen ineligible ID must produce exactly the predeclared structural
  abstention status and no factorized score row.
- An eligible-chain fallback, an ineligible-chain factorized success, a missing
  status, or any ID-set mismatch invalidates coverage.
- Paired factorized-versus-legacy resource measurements use the frozen eligible
  set; coverage is reported against all 663.

The runner may continue after an expected, frozen structural abstention. It must
still stop or invalidate the attempt for any unexpected fallback. Metrics remain
unreadable until this revised coverage validator succeeds.

### Reporting and promotion

The evaluator reports at least:

- `factorized_eligible_count / 663` and percentage;
- structural-abstention count and reason breakdown;
- factorized, legacy, Merizo, and Chainsaw metrics on the identical eligible set;
- full-set legacy and competitor metrics where their existing evidence permits;
- the nine factorized performance gates on the eligible set, clearly labeled as
  conditional rather than full-dataset estimates.

No abstained chain is silently dropped from coverage or replaced by its legacy
metric. Conditional accuracy and coverage are separate quantities.

The current model remains non-promotable as default-on if any of the following is
true:

- structural coverage is below 100%;
- any eligible chain falls back;
- any existing performance or resource gate fails;
- the previously frozen cache contract remains default-promotion incompatible.

An experimental opt-in may remain available because its structured fallback is
explicit. This design does not itself authorize default rollout.

### Audit disclosure

The earlier stopped CATH-663 attempt remains immutable audit evidence. Its
factorized failure status motivated this generic infrastructure correction, but
no CATH-663 metric was inspected. The final report must state that the protocol
was amended after a status-level context failure and must not describe the
replacement as an untouched first attempt.

## Failure semantics

Stable categories are:

- `structural_quality_abstention`: input residue evidence violates the strict
  complete-backbone policy;
- `feature_missing_context`: structurally eligible input reached the factorized
  pipeline but required typed evidence was unavailable;
- `identity`: candidate/lattice identity did not reconcile;
- existing candidate-local exclusions: ranking proceeded after invalid
  candidates were pruned transactionally.

Only the first category is allowed as an expected locked abstention, and only for
IDs frozen as ineligible before model execution. All other whole-chain fallbacks
remain invalid locked evidence.

## Verification strategy

### Rust unit and integration tests

- A complete standard residue is eligible.
- Missing N, CA, C, or O is reported in canonical backbone order.
- NaN or infinite required coordinates are ineligible.
- Blank and `A` alternate locations match DSSP semantics; other alternate
  locations alone do not satisfy eligibility.
- HETATM-only, nonstandard, and insertion-code residues do not enter coverage.
- Empty candidate-residue populations produce the stable empty reason.
- Inspector output is byte-identical across inspection and normal-run paths.
- Ineligible factorized input never reaches model inference and emits the exact
  structured abstention status.
- The same input under legacy selection preserves its existing selected result.
- Eligible fixtures still reach factorized selection.
- Missing Peeling evidence remains `feature_missing_context`, not structural
  abstention.

### Python and locked-protocol tests

- Input-only eligibility manifests reject duplicate IDs, reordered IDs, changed
  structure bytes, changed quality bytes, wrong policy, wrong binary, and
  eligible/ineligible overlap or omission.
- Coverage accepts only exact predeclared structural abstentions and exact actual
  factorized successes.
- Coverage rejects unexpected fallback, missing status, score row for an
  ineligible ID, absent score row for an eligible ID, and any hash mutation.
- Expected abstentions do not cause collection to stop; unexpected fallbacks do.
- Metric files are not parsed before coverage validation succeeds.
- Conditional comparisons use exactly the frozen eligible ID set for every tool.
- Coverage denominator remains 663 even when conditional metric denominators are
  smaller.

### Development-only validation

Before any replacement CATH-663 inspection:

- run focused Rust and Python suites;
- run `cargo check` and full `cargo test`;
- verify ordinary legacy behavior on representative complete and incomplete
  CATH-17287 inputs;
- verify known development examples: a complete chain succeeds, missing-CA and
  missing-O chains abstain, and empty-Peeling evidence retains its distinct code;
- regenerate the CATH-17287 input-only quality profile and reconcile its counts
  with the recorded diagnostic evidence;
- perform a fresh release build and freeze all source, binary, model, policy, and
  tool hashes.

No CATH-663 score, aggregate, competitor output, or model-performance metric may
be inspected during implementation or development validation.

## Future gap-aware model

A gap-tolerant successor requires a separate design and version. At minimum it
must define exact observed-residue set semantics, prevent endpoint ranges from
implicitly filling missing residues, include missingness information available at
inference time, train or calibrate with development-only synthetic deletion
patterns, predeclare acceptance/abstention thresholds, and receive a newly frozen
evaluation. Version 1 strict abstention must not be relaxed by a local heuristic.

## Completion criteria

This repair is complete when:

1. runtime and inspection paths share one tested structural-quality predicate;
2. factorized inference never runs on an ineligible chain;
3. abstention is machine-readable and distinct from factorized success;
4. legacy selection remains behaviorally unchanged;
5. the revised locked protocol preserves all 663 identities while comparing
   scores only on the frozen eligible subset;
6. all focused and full verification gates pass;
7. replacement intent and runtime artifacts are frozen before any further
   CATH-663 inspection;
8. the eventual report discloses coverage and the protocol amendment separately
   from conditional accuracy.
