use std::collections::BTreeMap;

use serde::Deserialize;

#[path = "../src/sword/factorized_ranker/generated_model.rs"]
mod generated_model;
#[path = "../src/sword/factorized_ranker/model.rs"]
mod model;
#[path = "../src/sword/factorized_ranker/partition.rs"]
mod partition;
#[path = "../src/sword/factorized_ranker/schema.rs"]
mod schema;

const GOLDEN: &str = include_str!("../../benchmark/models/factorized_ranker_v1_golden.json");
const TOLERANCE: f64 = 1.0e-10;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Golden {
    schema_version: u32,
    input_dtype: String,
    threshold_policy: String,
    retained_feature_families: Vec<String>,
    heads: Heads,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Heads {
    count: Head,
    candidate: Head,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Head {
    feature_names_sha256: String,
    representative: Representative,
    cases: Vec<Case>,
    errors: Vec<Mutation>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Representative {
    feature_index: usize,
    feature_name: String,
    following: f64,
    node_index: usize,
    previous: f64,
    raw_index: usize,
    raw_kind: String,
    threshold: f64,
    tree_index: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    kind: String,
    chain_id: String,
    shared_names: Vec<String>,
    item_names: Vec<String>,
    pair_feature_names: Vec<String>,
    shared_values: Vec<f64>,
    left_id: String,
    right_id: String,
    left_values: Vec<f64>,
    right_values: Vec<f64>,
    forward_pair_vector: Vec<f64>,
    reverse_pair_vector: Vec<f64>,
    forward_canonical_vector: Vec<f64>,
    reverse_canonical_vector: Vec<f64>,
    forward_tree_leaf_values: Vec<f64>,
    reverse_tree_leaf_values: Vec<f64>,
    forward_raw_margin: f64,
    reverse_raw_margin: f64,
    forward_probability: f64,
    reverse_probability: f64,
    symmetrized_probability: f64,
    borda_scores: BTreeMap<String, f64>,
    expected_equality_branch: String,
    tie_input: TieInput,
    winner: Winner,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", deny_unknown_fields)]
enum TieInput {
    #[serde(rename = "count")]
    Count {
        left_count: usize,
        right_count: usize,
        legacy_count: usize,
    },
    #[serde(rename = "candidate")]
    Candidate { tie_break: String },
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
#[serde(untagged)]
enum Winner {
    Count(usize),
    Candidate(String),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Mutation {
    kind: String,
    vector: Option<String>,
    index: Option<usize>,
    drop_last: Option<usize>,
}

#[derive(Clone, Copy)]
enum HeadKind {
    Count,
    Candidate,
}

fn assert_close(actual: f64, expected: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= TOLERANCE,
        "{label}: actual={actual:?}, expected={expected:?}, delta={:?}",
        (actual - expected).abs()
    );
}

fn assert_bits(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label} length");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{label}[{index}]: actual={actual:?}, expected={expected:?}"
        );
    }
}

// serde_json's default fast float parser can choose the adjacent f64 for a
// handful of halfway-looking decimal literals. The frozen artifact promises
// exact binary64 vectors, so parse those flat arrays from their original JSON
// lexemes with Rust's correctly rounded f64 parser before comparing bits.
fn exact_numeric_array(key: &str, occurrence: usize) -> Vec<f64> {
    let needle = format!("\"{key}\":[");
    let start = GOLDEN
        .match_indices(&needle)
        .nth(occurrence)
        .map(|(index, _)| index + needle.len())
        .unwrap_or_else(|| panic!("missing {key} occurrence {occurrence}"));
    let end = GOLDEN[start..]
        .find(']')
        .map(|offset| start + offset)
        .unwrap_or_else(|| panic!("unterminated {key} occurrence {occurrence}"));
    let body = &GOLDEN[start..end];
    if body.is_empty() {
        return Vec::new();
    }
    body.split(',')
        .map(|value| value.parse::<f64>().unwrap())
        .collect()
}

fn raw_case_ordinal(kind: HeadKind, case_index: usize) -> usize {
    // Canonical JSON sorts the two head keys: candidate's four cases precede
    // count's four cases in the immutable artifact.
    match kind {
        HeadKind::Candidate => case_index,
        HeadKind::Count => 4 + case_index,
    }
}

fn assemble(shared: &[f64], left: &[f64], right: &[f64]) -> Result<Vec<f64>, model::ModelError> {
    if left.len() != right.len() {
        return Err(model::ModelError::InputLength {
            expected: left.len(),
            actual: right.len(),
        });
    }
    let mut values = Vec::with_capacity(shared.len() + 2 * left.len());
    values.extend_from_slice(shared);
    for (left, right) in left.iter().zip(right) {
        values.push(left - right);
    }
    for (left, right) in left.iter().zip(right) {
        values.push((left - right).abs());
    }
    Ok(values)
}

fn model_for(kind: HeadKind) -> &'static model::StaticBoostedModel {
    match kind {
        HeadKind::Count => &generated_model::COUNT_MODEL,
        HeadKind::Candidate => &generated_model::CANDIDATE_MODEL,
    }
}

fn names_for(kind: HeadKind) -> &'static [&'static str] {
    match kind {
        HeadKind::Count => &generated_model::COUNT_FEATURE_NAMES,
        HeadKind::Candidate => &generated_model::CANDIDATE_FEATURE_NAMES,
    }
}

fn assert_case(kind: HeadKind, representative: &Representative, case: &Case, case_index: usize) {
    let model = model_for(kind);
    let generated_names = names_for(kind);
    let ordinal = raw_case_ordinal(kind, case_index);
    let shared_values = exact_numeric_array("shared_values", ordinal);
    let left_values = exact_numeric_array("left_values", ordinal);
    let right_values = exact_numeric_array("right_values", ordinal);
    let expected_forward_pair = exact_numeric_array("forward_pair_vector", ordinal);
    let expected_reverse_pair = exact_numeric_array("reverse_pair_vector", ordinal);
    let expected_forward_canonical = exact_numeric_array("forward_canonical_vector", ordinal);
    let expected_reverse_canonical = exact_numeric_array("reverse_canonical_vector", ordinal);
    let expected_forward_leaves = exact_numeric_array("forward_tree_leaf_values", ordinal);
    let expected_reverse_leaves = exact_numeric_array("reverse_tree_leaf_values", ordinal);
    assert_eq!(case.pair_feature_names.len(), generated_names.len());
    assert!(case
        .pair_feature_names
        .iter()
        .zip(generated_names)
        .all(|(golden, generated)| golden == generated));
    assert_eq!(case.shared_names.len(), shared_values.len());
    assert_eq!(case.item_names.len(), left_values.len());
    assert_eq!(case.item_names.len(), right_values.len());
    assert_eq!(case.shared_values.len(), shared_values.len());
    assert_eq!(case.left_values.len(), left_values.len());
    assert_eq!(case.right_values.len(), right_values.len());
    assert_eq!(case.forward_pair_vector.len(), expected_forward_pair.len());
    assert_eq!(case.reverse_pair_vector.len(), expected_reverse_pair.len());
    assert_eq!(
        case.forward_canonical_vector.len(),
        expected_forward_canonical.len()
    );
    assert_eq!(
        case.reverse_canonical_vector.len(),
        expected_reverse_canonical.len()
    );
    assert_eq!(
        case.forward_tree_leaf_values.len(),
        expected_forward_leaves.len()
    );
    assert_eq!(
        case.reverse_tree_leaf_values.len(),
        expected_reverse_leaves.len()
    );
    let expected_names = case
        .shared_names
        .iter()
        .cloned()
        .chain(case.item_names.iter().map(|name| format!("diff__{name}")))
        .chain(
            case.item_names
                .iter()
                .map(|name| format!("abs_diff__{name}")),
        )
        .collect::<Vec<_>>();
    assert_eq!(expected_names, case.pair_feature_names);

    let forward = assemble(&shared_values, &left_values, &right_values).unwrap();
    let reverse = assemble(&shared_values, &right_values, &left_values).unwrap();
    assert_bits(&forward, &expected_forward_pair, "forward pair");
    assert_bits(&reverse, &expected_reverse_pair, "reverse pair");

    let forward_canonical = model::canonicalize_input(&forward).unwrap();
    let reverse_canonical = model::canonicalize_input(&reverse).unwrap();
    assert_bits(
        &forward_canonical,
        &expected_forward_canonical,
        "forward canonical",
    );
    assert_bits(
        &reverse_canonical,
        &expected_reverse_canonical,
        "reverse canonical",
    );

    let forward_leaves = model::tree_leaf_values(model, &forward).unwrap();
    let reverse_leaves = model::tree_leaf_values(model, &reverse).unwrap();
    assert_bits(&forward_leaves, &expected_forward_leaves, "forward leaves");
    assert_bits(&reverse_leaves, &expected_reverse_leaves, "reverse leaves");

    let forward_raw = model::raw_score(model, &forward).unwrap();
    let reverse_raw = model::raw_score(model, &reverse).unwrap();
    assert_close(forward_raw, case.forward_raw_margin, "forward raw");
    assert_close(reverse_raw, case.reverse_raw_margin, "reverse raw");
    let forward_probability = model::predict_probability(model, &forward).unwrap();
    let reverse_probability = model::predict_probability(model, &reverse).unwrap();
    assert_close(
        forward_probability,
        case.forward_probability,
        "forward probability",
    );
    assert_close(
        reverse_probability,
        case.reverse_probability,
        "reverse probability",
    );
    let symmetrized = 0.5 * (forward_probability + 1.0 - reverse_probability);
    assert_close(
        symmetrized,
        case.symmetrized_probability,
        "symmetrized probability",
    );

    let identities = vec![case.right_id.clone(), case.left_id.clone()];
    let scores = model::normalized_borda(&identities, |left, right| {
        if left == &case.left_id && right == &case.right_id {
            model::predict_probability(model, &forward)
        } else if left == &case.right_id && right == &case.left_id {
            model::predict_probability(model, &reverse)
        } else {
            unreachable!("unexpected golden identity")
        }
    })
    .unwrap()
    .into_iter()
    .collect::<BTreeMap<_, _>>();
    assert_close(scores[&case.left_id], symmetrized, "left Borda");
    assert_close(scores[&case.right_id], 1.0 - symmetrized, "right Borda");

    let (left_key, right_key, tie_winner) = match &case.tie_input {
        TieInput::Count {
            left_count,
            right_count,
            legacy_count,
        } => {
            let exact_tie = if legacy_count == left_count || legacy_count == right_count {
                *legacy_count
            } else {
                (*left_count).min(*right_count)
            };
            (
                left_count.to_string(),
                right_count.to_string(),
                Winner::Count(exact_tie),
            )
        }
        TieInput::Candidate { tie_break } => {
            assert_eq!(tie_break, "lexical");
            (
                case.left_id.clone(),
                case.right_id.clone(),
                Winner::Candidate(std::cmp::min(&case.left_id, &case.right_id).clone()),
            )
        }
    };
    assert_close(
        case.borda_scores[&left_key],
        symmetrized,
        "golden left Borda",
    );
    assert_close(
        case.borda_scores[&right_key],
        1.0 - symmetrized,
        "golden right Borda",
    );
    let actual_winner = if symmetrized > 0.5 {
        match &case.tie_input {
            TieInput::Count { left_count, .. } => Winner::Count(*left_count),
            TieInput::Candidate { .. } => Winner::Candidate(case.left_id.clone()),
        }
    } else if symmetrized < 0.5 {
        match &case.tie_input {
            TieInput::Count { right_count, .. } => Winner::Count(*right_count),
            TieInput::Candidate { .. } => Winner::Candidate(case.right_id.clone()),
        }
    } else {
        tie_winner
    };
    assert_eq!(actual_winner, case.winner);

    if case.kind == "threshold_equal" {
        assert_eq!(case.expected_equality_branch, "left");
        assert_eq!(
            representative.feature_name,
            generated_names[representative.feature_index]
        );
        assert_eq!(representative.tree_index, 0);
        assert_eq!(representative.node_index, 0);
        assert_eq!(representative.raw_kind, "signed_difference");
        assert_eq!(
            usize::from(model.nodes[representative.node_index].feature),
            representative.feature_index
        );
        assert_eq!(
            model.nodes[representative.node_index].threshold.to_bits(),
            representative.threshold.to_bits()
        );
        assert_eq!(
            usize::from(model.trees[representative.tree_index].root),
            representative.node_index
        );
        assert_eq!(
            expected_forward_canonical[representative.feature_index].to_bits(),
            representative.threshold.to_bits()
        );
        assert!(
            expected_forward_canonical[representative.feature_index]
                <= model.nodes[representative.node_index].threshold
        );
        assert_eq!(model.nodes[representative.node_index].left, 1);
    } else {
        assert_eq!(case.expected_equality_branch, "not_applicable");
    }
    match case.kind.as_str() {
        "threshold_below" => assert_eq!(
            forward[representative.feature_index].to_bits(),
            representative.previous.to_bits()
        ),
        "threshold_equal" => assert_eq!(
            forward[representative.feature_index].to_bits(),
            representative.threshold.to_bits()
        ),
        "threshold_above" => assert_eq!(
            forward[representative.feature_index].to_bits(),
            representative.following.to_bits()
        ),
        "realistic" => {}
        _ => panic!("unexpected golden case kind"),
    }
}

fn assert_mutations(kind: HeadKind, head: &Head) {
    let model = model_for(kind);
    let case = &head.cases[0];
    let expected_kinds = [
        "nan_at",
        "positive_infinity_at",
        "negative_infinity_at",
        "length_mismatch",
        "duplicate_item_id",
        "float32_overflow_at",
    ];
    assert_eq!(head.errors.len(), expected_kinds.len());
    for (mutation, expected_kind) in head.errors.iter().zip(expected_kinds) {
        assert_eq!(mutation.kind, expected_kind);
        match mutation.kind.as_str() {
            "nan_at" | "positive_infinity_at" | "negative_infinity_at" | "float32_overflow_at" => {
                let index = mutation.index.unwrap();
                let vector = mutation.vector.as_deref().unwrap();
                let mut shared = case.shared_values.clone();
                let mut left = case.left_values.clone();
                let mut right = case.right_values.clone();
                let replacement = match mutation.kind.as_str() {
                    "nan_at" => f64::NAN,
                    "positive_infinity_at" => f64::INFINITY,
                    "negative_infinity_at" => f64::NEG_INFINITY,
                    "float32_overflow_at" => f64::MAX,
                    _ => unreachable!(),
                };
                match vector {
                    "shared" => shared[index] = replacement,
                    "left" => left[index] = replacement,
                    "right" => right[index] = replacement,
                    _ => panic!("unknown mutation vector"),
                }
                let pair = assemble(&shared, &left, &right).unwrap();
                let error = model::predict_probability(model, &pair).unwrap_err();
                match mutation.kind.as_str() {
                    "float32_overflow_at" => {
                        assert!(matches!(error, model::ModelError::Float32Overflow { .. }))
                    }
                    _ => assert!(matches!(error, model::ModelError::NonFiniteInput { .. })),
                }
            }
            "length_mismatch" => {
                assert_eq!(mutation.vector.as_deref(), Some("right"));
                let mut pair = case.forward_pair_vector.clone();
                for _ in 0..mutation.drop_last.unwrap() {
                    pair.pop();
                }
                assert!(matches!(
                    model::predict_probability(model, &pair),
                    Err(model::ModelError::InputLength { .. })
                ));
            }
            "duplicate_item_id" => {
                assert!(matches!(
                    model::normalized_borda(
                        &[case.left_id.clone(), case.left_id.clone()],
                        |_, _| { Ok(0.5) }
                    ),
                    Err(model::ModelError::DuplicateIdentity { .. })
                ));
            }
            _ => unreachable!(),
        }
    }
}

fn assert_head(kind: HeadKind, head: &Head) {
    assert_eq!(head.feature_names_sha256.len(), 64);
    assert_eq!(head.cases.len(), 4);
    assert_eq!(
        head.cases
            .iter()
            .map(|case| case.kind.as_str())
            .collect::<Vec<_>>(),
        [
            "realistic",
            "threshold_below",
            "threshold_equal",
            "threshold_above"
        ]
    );
    for (case_index, case) in head.cases.iter().enumerate() {
        assert!(!case.chain_id.is_empty());
        assert_case(kind, &head.representative, case, case_index);
    }
    assert!(head.representative.previous < head.representative.threshold);
    assert!(head.representative.following > head.representative.threshold);
    assert!(head.representative.raw_index < head.cases[0].item_names.len());
    assert_mutations(kind, head);
}

#[test]
fn frozen_python_golden_matches_rust_tree_runtime() {
    let golden: Golden = serde_json::from_str(GOLDEN).unwrap();
    assert_eq!(golden.schema_version, schema::FEATURE_SCHEMA_VERSION);
    assert_eq!(golden.input_dtype, generated_model::MODEL_INPUT_DTYPE);
    assert_eq!(
        golden.threshold_policy,
        generated_model::MODEL_THRESHOLD_POLICY
    );
    assert_eq!(
        golden.retained_feature_families,
        generated_model::RETAINED_FEATURE_FAMILIES
    );
    model::validate_model(&generated_model::COUNT_MODEL).unwrap();
    model::validate_model(&generated_model::CANDIDATE_MODEL).unwrap();
    assert_head(HeadKind::Count, &golden.heads.count);
    assert_head(HeadKind::Candidate, &golden.heads.candidate);
}
