use std::collections::BTreeSet;

use super::schema::FEATURE_SCHEMA_VERSION;

pub(crate) const LEAF_FEATURE: u16 = u16::MAX;

#[derive(Debug, Clone, Copy)]
pub(crate) struct StaticNode {
    pub feature: u16,
    pub threshold: f64,
    pub left: u16,
    pub right: u16,
    pub leaf_value: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct StaticTree {
    pub root: u16,
}

pub(crate) struct StaticBoostedModel {
    pub schema_version: u32,
    pub feature_names: &'static [&'static str],
    pub initial_log_odds: f64,
    pub learning_rate: f64,
    pub trees: &'static [StaticTree],
    pub nodes: &'static [StaticNode],
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum ModelError {
    #[error("model schema version {actual} does not match expected {expected}")]
    SchemaVersion { expected: u32, actual: u32 },
    #[error("model policy {policy} is {actual}, expected {expected}")]
    Policy {
        policy: &'static str,
        expected: &'static str,
        actual: &'static str,
    },
    #[error("invalid {head} model feature schema")]
    FeatureSchema { head: &'static str },
    #[error("empty model feature name at index {index}")]
    EmptyFeatureName { index: usize },
    #[error("duplicate model feature name at index {index}")]
    DuplicateFeatureName { index: usize },
    #[error("model forest is empty")]
    EmptyForest,
    #[error("model has {actual} trees, maximum is {maximum}")]
    TreeLimit { actual: usize, maximum: usize },
    #[error("models have {actual} nodes, maximum is {maximum}")]
    NodeLimit { actual: usize, maximum: usize },
    #[error("nonfinite model field {field}")]
    NonFiniteModel { field: &'static str },
    #[error("model learning rate must be positive")]
    NonPositiveLearningRate,
    #[error("tree {tree} root {root} is invalid")]
    InvalidRoot { tree: usize, root: u16 },
    #[error("node {node} has nonfinite field {field}")]
    NonFiniteNode { node: usize, field: &'static str },
    #[error("leaf node {node} has invalid sentinels")]
    InvalidLeafSentinel { node: usize },
    #[error("internal node {node} has invalid sentinels")]
    InvalidInternalSentinel { node: usize },
    #[error("internal node {node} feature {feature} is invalid")]
    InvalidFeature { node: usize, feature: u16 },
    #[error("node {node} child {child} is invalid")]
    InvalidChild { node: usize, child: u16 },
    #[error("node {node} uses child {child} twice")]
    DuplicateChild { node: usize, child: u16 },
    #[error("node {node} threshold is not canonical float32")]
    NonCanonicalThreshold { node: usize },
    #[error("node {node} links to itself")]
    SelfLink { node: usize },
    #[error("cycle reaches node {node}")]
    Cycle { node: usize },
    #[error("node {node} is shared")]
    SharedNode { node: usize },
    #[error("tree {tree} reuses root {root}")]
    ReusedRoot { tree: usize, root: u16 },
    #[error("tree {tree} reaches node {node} at depth {depth}")]
    DepthExceeded {
        tree: usize,
        node: usize,
        depth: usize,
    },
    #[error("node {node} is unreachable")]
    UnreachableNode { node: usize },
    #[error("model input length {actual} does not match expected {expected}")]
    InputLength { expected: usize, actual: usize },
    #[error("model input {index} is nonfinite")]
    NonFiniteInput { index: usize },
    #[error("model input {index} overflows float32")]
    Float32Overflow { index: usize },
    #[error("model output {stage} is nonfinite")]
    NonFiniteOutput { stage: &'static str },
    #[error("invalid probability for pair ({left}, {right})")]
    InvalidProbability { left: usize, right: usize },
    #[error("ranking group is empty")]
    EmptyGroup,
    #[error("duplicate identity at sorted index {index}")]
    DuplicateIdentity { index: usize },
    #[error("factorized lattice group mismatch")]
    LatticeGroupMismatch,
    #[error("factorized lattice identity mismatch")]
    LatticeIdentityMismatch,
}

pub(crate) fn validate_model(model: &StaticBoostedModel) -> Result<(), ModelError> {
    if model.schema_version != FEATURE_SCHEMA_VERSION {
        return Err(ModelError::SchemaVersion {
            expected: FEATURE_SCHEMA_VERSION,
            actual: model.schema_version,
        });
    }
    if model.feature_names.is_empty() {
        return Err(ModelError::FeatureSchema { head: "model" });
    }
    let mut feature_names = BTreeSet::new();
    for (index, name) in model.feature_names.iter().enumerate() {
        if name.is_empty() {
            return Err(ModelError::EmptyFeatureName { index });
        }
        if !feature_names.insert(*name) {
            return Err(ModelError::DuplicateFeatureName { index });
        }
    }
    if model.trees.is_empty() {
        return Err(ModelError::EmptyForest);
    }
    if model.trees.len() > 96 {
        return Err(ModelError::TreeLimit {
            actual: model.trees.len(),
            maximum: 96,
        });
    }
    if !model.initial_log_odds.is_finite() {
        return Err(ModelError::NonFiniteModel {
            field: "initial_log_odds",
        });
    }
    if !model.learning_rate.is_finite() {
        return Err(ModelError::NonFiniteModel {
            field: "learning_rate",
        });
    }
    if model.learning_rate <= 0.0 {
        return Err(ModelError::NonPositiveLearningRate);
    }

    for (index, node) in model.nodes.iter().enumerate() {
        if !node.threshold.is_finite() {
            return Err(ModelError::NonFiniteNode {
                node: index,
                field: "threshold",
            });
        }
        if !node.leaf_value.is_finite() {
            return Err(ModelError::NonFiniteNode {
                node: index,
                field: "leaf_value",
            });
        }
        if node.feature == LEAF_FEATURE {
            if node.left != LEAF_FEATURE
                || node.right != LEAF_FEATURE
                || node.threshold.to_bits() != 0.0_f64.to_bits()
            {
                return Err(ModelError::InvalidLeafSentinel { node: index });
            }
        } else {
            if node.left == LEAF_FEATURE
                || node.right == LEAF_FEATURE
                || node.leaf_value.to_bits() != 0.0_f64.to_bits()
            {
                return Err(ModelError::InvalidInternalSentinel { node: index });
            }
            if usize::from(node.feature) >= model.feature_names.len() {
                return Err(ModelError::InvalidFeature {
                    node: index,
                    feature: node.feature,
                });
            }
            for child in [node.left, node.right] {
                if usize::from(child) >= model.nodes.len() {
                    return Err(ModelError::InvalidChild { node: index, child });
                }
                if usize::from(child) == index {
                    return Err(ModelError::SelfLink { node: index });
                }
            }
            if node.left == node.right {
                return Err(ModelError::DuplicateChild {
                    node: index,
                    child: node.left,
                });
            }
            if ((node.threshold as f32) as f64).to_bits() != node.threshold.to_bits() {
                return Err(ModelError::NonCanonicalThreshold { node: index });
            }
        }
    }

    fn visit(
        model: &StaticBoostedModel,
        tree: usize,
        node_index: usize,
        depth: usize,
        colors: &mut [u8],
        owners: &mut [Option<usize>],
    ) -> Result<(), ModelError> {
        if depth > 3 {
            return Err(ModelError::DepthExceeded {
                tree,
                node: node_index,
                depth,
            });
        }
        if colors[node_index] == 1 {
            return Err(ModelError::Cycle { node: node_index });
        }
        if colors[node_index] == 2 || owners[node_index].is_some() {
            return Err(ModelError::SharedNode { node: node_index });
        }
        colors[node_index] = 1;
        owners[node_index] = Some(tree);
        let node = model.nodes[node_index];
        if node.feature != LEAF_FEATURE {
            visit(
                model,
                tree,
                usize::from(node.left),
                depth + 1,
                colors,
                owners,
            )?;
            visit(
                model,
                tree,
                usize::from(node.right),
                depth + 1,
                colors,
                owners,
            )?;
        }
        colors[node_index] = 2;
        Ok(())
    }

    let mut colors = vec![0_u8; model.nodes.len()];
    let mut owners = vec![None; model.nodes.len()];
    for (tree, tree_data) in model.trees.iter().enumerate() {
        let root = usize::from(tree_data.root);
        if root >= model.nodes.len() {
            return Err(ModelError::InvalidRoot {
                tree,
                root: tree_data.root,
            });
        }
        if owners[root].is_some() {
            return Err(ModelError::ReusedRoot {
                tree,
                root: tree_data.root,
            });
        }
        visit(model, tree, root, 0, &mut colors, &mut owners)?;
    }
    if let Some(node) = owners.iter().position(Option::is_none) {
        return Err(ModelError::UnreachableNode { node });
    }
    Ok(())
}

pub(crate) fn canonicalize_input(features: &[f64]) -> Result<Vec<f64>, ModelError> {
    features
        .iter()
        .enumerate()
        .map(|(index, value)| {
            if !value.is_finite() {
                return Err(ModelError::NonFiniteInput { index });
            }
            let narrowed = *value as f32;
            if !narrowed.is_finite() {
                return Err(ModelError::Float32Overflow { index });
            }
            Ok(narrowed as f64)
        })
        .collect()
}

pub(crate) fn tree_leaf_values(
    model: &StaticBoostedModel,
    features: &[f64],
) -> Result<Vec<f64>, ModelError> {
    validate_model(model)?;
    let features = prepare_input(model, features)?;
    Ok(tree_leaf_values_validated(model, &features))
}

pub(crate) fn raw_score(model: &StaticBoostedModel, features: &[f64]) -> Result<f64, ModelError> {
    validate_model(model)?;
    let features = prepare_input(model, features)?;
    raw_score_validated(model, &features)
}

pub(crate) fn predict_probability(
    model: &StaticBoostedModel,
    features: &[f64],
) -> Result<f64, ModelError> {
    validate_model(model)?;
    predict_probability_validated(model, features)
}

pub(crate) fn normalized_borda<I, F>(
    identities: &[I],
    mut directional_probability: F,
) -> Result<Vec<(I, f64)>, ModelError>
where
    I: Ord + Clone,
    F: FnMut(&I, &I) -> Result<f64, ModelError>,
{
    if identities.is_empty() {
        return Err(ModelError::EmptyGroup);
    }
    let mut identities = identities.to_vec();
    identities.sort();
    if let Some(index) = identities.windows(2).position(|pair| pair[0] == pair[1]) {
        return Err(ModelError::DuplicateIdentity { index: index + 1 });
    }
    if identities.len() == 1 {
        return Ok(vec![(identities.remove(0), 1.0)]);
    }

    let mut totals = vec![0.0; identities.len()];
    for left in 0..identities.len() {
        for right in (left + 1)..identities.len() {
            let forward = directional_probability(&identities[left], &identities[right])?;
            validate_probability(forward, left, right)?;
            let reverse = directional_probability(&identities[right], &identities[left])?;
            validate_probability(reverse, right, left)?;
            let left_probability = 0.5 * (forward + 1.0 - reverse);
            validate_probability(left_probability, left, right)?;
            let right_probability = 1.0 - left_probability;
            validate_probability(right_probability, right, left)?;
            totals[left] += left_probability;
            totals[right] += right_probability;
        }
    }
    let denominator = (identities.len() - 1) as f64;
    identities
        .into_iter()
        .zip(totals)
        .enumerate()
        .map(|(index, (identity, total))| {
            let score = total / denominator;
            validate_probability(score, index, index)?;
            Ok((identity, score))
        })
        .collect()
}

fn prepare_input(model: &StaticBoostedModel, features: &[f64]) -> Result<Vec<f64>, ModelError> {
    if features.len() != model.feature_names.len() {
        return Err(ModelError::InputLength {
            expected: model.feature_names.len(),
            actual: features.len(),
        });
    }
    canonicalize_input(features)
}

fn tree_leaf_values_validated(model: &StaticBoostedModel, features: &[f64]) -> Vec<f64> {
    let mut leaves = Vec::with_capacity(model.trees.len());
    for tree in model.trees {
        let mut node_index = usize::from(tree.root);
        loop {
            let node = model.nodes[node_index];
            if node.feature == LEAF_FEATURE {
                leaves.push(node.leaf_value);
                break;
            }
            node_index = usize::from(if features[usize::from(node.feature)] <= node.threshold {
                node.left
            } else {
                node.right
            });
        }
    }
    leaves
}

fn raw_score_validated(model: &StaticBoostedModel, features: &[f64]) -> Result<f64, ModelError> {
    let mut leaf_sum = 0.0;
    for leaf in tree_leaf_values_validated(model, features) {
        leaf_sum += leaf;
    }
    let raw = model.initial_log_odds + model.learning_rate * leaf_sum;
    if !raw.is_finite() {
        return Err(ModelError::NonFiniteOutput {
            stage: "raw_margin",
        });
    }
    Ok(raw)
}

pub(crate) fn predict_probability_validated(
    model: &StaticBoostedModel,
    features: &[f64],
) -> Result<f64, ModelError> {
    let features = prepare_input(model, features)?;
    let raw = raw_score_validated(model, &features)?;
    let probability = if raw >= 0.0 {
        1.0 / (1.0 + (-raw).exp())
    } else {
        let exp_raw = raw.exp();
        exp_raw / (1.0 + exp_raw)
    };
    if !probability.is_finite() {
        return Err(ModelError::NonFiniteOutput {
            stage: "probability",
        });
    }
    validate_probability(probability, 0, 0)?;
    Ok(probability)
}

fn validate_probability(probability: f64, left: usize, right: usize) -> Result<(), ModelError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(ModelError::InvalidProbability { left, right });
    }
    Ok(())
}

#[cfg(test)]
mod factorized_model_tests {
    use super::*;
    use std::cell::RefCell;

    static FEATURES: &[&str] = &["x"];

    fn one_split_model(threshold: f64, left_value: f64, right_value: f64) -> StaticBoostedModel {
        let nodes = Box::leak(Box::new([
            StaticNode {
                feature: 0,
                threshold,
                left: 1,
                right: 2,
                leaf_value: 0.0,
            },
            StaticNode {
                feature: LEAF_FEATURE,
                threshold: 0.0,
                left: LEAF_FEATURE,
                right: LEAF_FEATURE,
                leaf_value: left_value,
            },
            StaticNode {
                feature: LEAF_FEATURE,
                threshold: 0.0,
                left: LEAF_FEATURE,
                right: LEAF_FEATURE,
                leaf_value: right_value,
            },
        ]));
        let trees = Box::leak(Box::new([StaticTree { root: 0 }]));
        StaticBoostedModel {
            schema_version: super::super::schema::FEATURE_SCHEMA_VERSION,
            feature_names: FEATURES,
            initial_log_odds: 0.0,
            learning_rate: 1.0,
            trees,
            nodes,
        }
    }

    #[test]
    fn equality_uses_the_left_child_after_float32_canonicalization() {
        let model = one_split_model(1.25, -1.0, 1.0);
        let equal = 1.25_f32;
        let previous = f32::from_bits(equal.to_bits() - 1) as f64;
        let following = f32::from_bits(equal.to_bits() + 1) as f64;
        assert_eq!(tree_leaf_values(&model, &[previous]).unwrap(), vec![-1.0]);
        assert_eq!(
            tree_leaf_values(&model, &[equal as f64]).unwrap(),
            vec![-1.0]
        );
        assert_eq!(tree_leaf_values(&model, &[following]).unwrap(), vec![1.0]);
    }

    #[test]
    fn finite_f64_that_overflows_float32_is_rejected() {
        let model = one_split_model(0.0, -1.0, 1.0);
        assert_eq!(
            predict_probability(&model, &[f64::MAX]).unwrap_err(),
            ModelError::Float32Overflow { index: 0 }
        );
    }

    #[test]
    fn canonicalization_preserves_negative_zero_and_rejects_nonfinite_first() {
        let values = canonicalize_input(&[-0.0, 1.0 + f64::EPSILON]).unwrap();
        assert_eq!(values[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(values[1].to_bits(), 1.0_f64.to_bits());
        assert_eq!(
            canonicalize_input(&[f64::NAN]).unwrap_err(),
            ModelError::NonFiniteInput { index: 0 }
        );
    }

    #[test]
    fn ordered_tree_sum_and_stable_sigmoid_match_hand_calculation() {
        let mut model = one_split_model(0.0, 2.0, -2.0);
        model.initial_log_odds = -0.5;
        model.learning_rate = 0.25;
        assert_eq!(raw_score(&model, &[-0.0]).unwrap(), 0.0);
        assert_eq!(predict_probability(&model, &[-0.0]).unwrap(), 0.5);

        model.initial_log_odds = -1000.0;
        assert_eq!(predict_probability(&model, &[1.0]).unwrap(), 0.0);
        model.initial_log_odds = 1000.0;
        assert_eq!(predict_probability(&model, &[-1.0]).unwrap(), 1.0);

        let nodes = Box::leak(Box::new([
            StaticNode {
                feature: LEAF_FEATURE,
                threshold: 0.0,
                left: LEAF_FEATURE,
                right: LEAF_FEATURE,
                leaf_value: 1.0e16,
            },
            StaticNode {
                feature: LEAF_FEATURE,
                threshold: 0.0,
                left: LEAF_FEATURE,
                right: LEAF_FEATURE,
                leaf_value: -1.0e16,
            },
            StaticNode {
                feature: LEAF_FEATURE,
                threshold: 0.0,
                left: LEAF_FEATURE,
                right: LEAF_FEATURE,
                leaf_value: 1.0,
            },
        ]));
        let trees = Box::leak(Box::new([
            StaticTree { root: 0 },
            StaticTree { root: 1 },
            StaticTree { root: 2 },
        ]));
        let ordered = StaticBoostedModel {
            schema_version: super::super::schema::FEATURE_SCHEMA_VERSION,
            feature_names: FEATURES,
            initial_log_odds: 0.5,
            learning_rate: 0.25,
            trees,
            nodes,
        };
        assert_eq!(raw_score(&ordered, &[0.0]).unwrap(), 0.75);
    }

    #[test]
    fn borda_is_permutation_invariant_and_calls_each_orientation_once() {
        fn rank(order: &[usize]) -> (Vec<(usize, f64)>, Vec<(usize, usize)>) {
            let calls = RefCell::new(Vec::new());
            let ranked = normalized_borda(order, |left, right| {
                calls.borrow_mut().push((*left, *right));
                Ok(if left < right { 0.75 } else { 0.25 })
            })
            .unwrap();
            (ranked, calls.into_inner())
        }
        let (first, calls) = rank(&[2, 0, 1]);
        let (second, _) = rank(&[1, 2, 0]);
        assert_eq!(first, second);
        assert_eq!(calls, vec![(0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)]);
    }

    #[test]
    fn borda_handles_singletons_and_rejects_invalid_groups() {
        let calls = RefCell::new(0);
        let singleton = normalized_borda(&[7], |_, _| {
            *calls.borrow_mut() += 1;
            Ok(0.5)
        })
        .unwrap();
        assert_eq!(singleton, vec![(7, 1.0)]);
        assert_eq!(*calls.borrow(), 0);
        assert_eq!(
            normalized_borda::<usize, _>(&[], |_, _| Ok(0.5)).unwrap_err(),
            ModelError::EmptyGroup
        );
        assert_eq!(
            normalized_borda(&[2, 1, 2], |_, _| Ok(0.5)).unwrap_err(),
            ModelError::DuplicateIdentity { index: 2 }
        );
        assert!(matches!(
            normalized_borda(&[1, 2], |_, _| Ok(f64::NAN)),
            Err(ModelError::InvalidProbability { .. })
        ));
        assert!(matches!(
            normalized_borda(&[1, 2], |_, _| Ok(1.5)),
            Err(ModelError::InvalidProbability { .. })
        ));
    }

    #[test]
    fn malformed_forests_have_stable_errors() {
        let model = one_split_model(0.1, -1.0, 1.0);
        assert_eq!(
            validate_model(&model).unwrap_err(),
            ModelError::NonCanonicalThreshold { node: 0 }
        );

        let mut empty = one_split_model(0.0, -1.0, 1.0);
        empty.trees = &[];
        assert_eq!(validate_model(&empty).unwrap_err(), ModelError::EmptyForest);

        let mut empty_features = one_split_model(0.0, -1.0, 1.0);
        empty_features.feature_names = &[];
        assert_eq!(
            validate_model(&empty_features).unwrap_err(),
            ModelError::FeatureSchema { head: "model" }
        );

        let mut nonfinite = one_split_model(0.0, -1.0, 1.0);
        nonfinite.initial_log_odds = f64::NAN;
        assert_eq!(
            validate_model(&nonfinite).unwrap_err(),
            ModelError::NonFiniteModel {
                field: "initial_log_odds"
            }
        );

        let nonfinite_node = static_model(
            vec![
                StaticNode {
                    threshold: f64::INFINITY,
                    ..internal(1, 2)
                },
                leaf(0.0),
                leaf(0.0),
            ],
            vec![0],
        );
        assert_eq!(
            validate_model(&nonfinite_node).unwrap_err(),
            ModelError::NonFiniteNode {
                node: 0,
                field: "threshold"
            }
        );

        let invalid_leaf = Box::leak(Box::new([StaticNode {
            feature: LEAF_FEATURE,
            threshold: -0.0,
            left: LEAF_FEATURE,
            right: LEAF_FEATURE,
            leaf_value: 0.0,
        }]));
        let tree = Box::leak(Box::new([StaticTree { root: 0 }]));
        let invalid_leaf = StaticBoostedModel {
            schema_version: super::super::schema::FEATURE_SCHEMA_VERSION,
            feature_names: FEATURES,
            initial_log_odds: 0.0,
            learning_rate: 1.0,
            trees: tree,
            nodes: invalid_leaf,
        };
        assert_eq!(
            validate_model(&invalid_leaf).unwrap_err(),
            ModelError::InvalidLeafSentinel { node: 0 }
        );
    }

    fn leaf(value: f64) -> StaticNode {
        StaticNode {
            feature: LEAF_FEATURE,
            threshold: 0.0,
            left: LEAF_FEATURE,
            right: LEAF_FEATURE,
            leaf_value: value,
        }
    }

    fn internal(left: u16, right: u16) -> StaticNode {
        StaticNode {
            feature: 0,
            threshold: 0.0,
            left,
            right,
            leaf_value: 0.0,
        }
    }

    fn static_model(nodes: Vec<StaticNode>, roots: Vec<u16>) -> StaticBoostedModel {
        let nodes = Box::leak(nodes.into_boxed_slice());
        let trees = Box::leak(
            roots
                .into_iter()
                .map(|root| StaticTree { root })
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        );
        StaticBoostedModel {
            schema_version: super::super::schema::FEATURE_SCHEMA_VERSION,
            feature_names: FEATURES,
            initial_log_odds: 0.0,
            learning_rate: 1.0,
            trees,
            nodes,
        }
    }

    #[test]
    fn graph_validation_rejects_cycles_sharing_depth_and_unreachable_nodes() {
        let cycle = static_model(
            vec![internal(1, 2), internal(0, 3), leaf(0.0), leaf(0.0)],
            vec![0],
        );
        assert_eq!(
            validate_model(&cycle).unwrap_err(),
            ModelError::Cycle { node: 0 }
        );

        let shared = static_model(
            vec![
                internal(1, 2),
                internal(3, 4),
                internal(3, 5),
                leaf(0.0),
                leaf(0.0),
                leaf(0.0),
            ],
            vec![0],
        );
        assert_eq!(
            validate_model(&shared).unwrap_err(),
            ModelError::SharedNode { node: 3 }
        );

        let unreachable = static_model(
            vec![internal(1, 2), leaf(0.0), leaf(0.0), leaf(0.0)],
            vec![0],
        );
        assert_eq!(
            validate_model(&unreachable).unwrap_err(),
            ModelError::UnreachableNode { node: 3 }
        );

        let reused_root = static_model(vec![leaf(0.0), leaf(0.0)], vec![0, 0]);
        assert_eq!(
            validate_model(&reused_root).unwrap_err(),
            ModelError::ReusedRoot { tree: 1, root: 0 }
        );

        let too_deep = static_model(
            vec![
                internal(1, 5),
                internal(2, 6),
                internal(3, 7),
                internal(4, 8),
                leaf(0.0),
                leaf(0.0),
                leaf(0.0),
                leaf(0.0),
                leaf(0.0),
            ],
            vec![0],
        );
        assert_eq!(
            validate_model(&too_deep).unwrap_err(),
            ModelError::DepthExceeded {
                tree: 0,
                node: 4,
                depth: 4
            }
        );
    }

    #[test]
    fn overflowing_ordered_leaf_sum_is_a_nonfinite_output_error() {
        let model = static_model(vec![leaf(f64::MAX), leaf(f64::MAX)], vec![0, 1]);
        assert_eq!(
            raw_score(&model, &[0.0]).unwrap_err(),
            ModelError::NonFiniteOutput {
                stage: "raw_margin"
            }
        );
    }
}
