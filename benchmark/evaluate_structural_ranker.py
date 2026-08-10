"""Train and evaluate a standalone SWORD structural candidate ranker.

This evaluator uses only SWORD-produced candidate features: contact ratios,
domain-shape statistics, and segmentation topology.  It deliberately excludes
energy Z-scores and all external predictors.  The training table must come from
CATH-17287; CATH-663 is consumed only once as a held-out evaluation table.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import GroupShuffleSplit

from benchmark.candidate_geometry import CANDIDATE_GEOMETRY_FIELDS


FAST_BASE_FEATURES = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "modal_count_distance",
]
FEATURES = [*FAST_BASE_FEATURES, *CANDIDATE_GEOMETRY_FIELDS]


@dataclass(frozen=True)
class SelectionMetrics:
    n_chains: int
    ndo: float
    d_count_acc: float
    boundary_f1_10: float
    matched_dice: float
    predicted_domains: float


def load_table(path: Path) -> pd.DataFrame:
    """Load only candidates with complete ranker features and CATH labels."""
    frame = pd.read_csv(path)
    required = ["chain_id", "ndo", "d_count_acc", "boundary_f1_10", "matched_dice", *FEATURES]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    for column in ["ndo", "d_count_acc", "boundary_f1_10", "matched_dice", *FEATURES]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame.dropna(subset=required).reset_index(drop=True)


def make_model(seed: int) -> ExtraTreesRegressor:
    """A compact, regularized non-linear ranker with no runtime dependency."""
    return ExtraTreesRegressor(
        n_estimators=600,
        min_samples_leaf=4,
        max_features=0.8,
        n_jobs=-1,
        random_state=seed,
    )


def select_top1(frame: pd.DataFrame, predictions: np.ndarray) -> pd.DataFrame:
    """Pick the highest predicted candidate separately for each protein chain."""
    scored = frame.copy()
    scored["ranker_score"] = predictions
    indices = scored.groupby("chain_id", sort=False)["ranker_score"].idxmax()
    return scored.loc[indices].reset_index(drop=True)


def metrics(selected: pd.DataFrame) -> SelectionMetrics:
    return SelectionMetrics(
        n_chains=int(selected["chain_id"].nunique()),
        ndo=float(selected["ndo"].mean()),
        d_count_acc=float(selected["d_count_acc"].mean()),
        boundary_f1_10=float(selected["boundary_f1_10"].mean()),
        matched_dice=float(selected["matched_dice"].mean()),
        predicted_domains=float(selected["n_pred_domains"].mean()),
    )


def format_metrics(label: str, value: SelectionMetrics) -> str:
    return (
        f"{label}: chains={value.n_chains} NDO={value.ndo:.4f} "
        f"d_count_acc={value.d_count_acc:.4f} boundary_f1_10={value.boundary_f1_10:.4f} "
        f"matched_dice={value.matched_dice:.4f} predicted_domains={value.predicted_domains:.3f}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--selection-out", type=Path, default=None)
    parser.add_argument("--report-out", type=Path, default=None)
    args = parser.parse_args()

    train = load_table(args.train)
    test = load_table(args.test)
    train_x = train[FEATURES].to_numpy(dtype=float)
    train_y = train["ndo"].to_numpy(dtype=float)
    groups = train["chain_id"].to_numpy()

    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=args.validation_fraction,
        random_state=args.seed,
    )
    fit_index, validation_index = next(splitter.split(train_x, train_y, groups))
    validation_model = make_model(args.seed)
    validation_model.fit(train_x[fit_index], train_y[fit_index])
    validation_selected = select_top1(
        train.iloc[validation_index],
        validation_model.predict(train_x[validation_index]),
    )

    model = make_model(args.seed)
    model.fit(train_x, train_y)
    test_selected = select_top1(
        test,
        model.predict(test[FEATURES].to_numpy(dtype=float)),
    )

    validation_metrics = metrics(validation_selected)
    test_metrics = metrics(test_selected)
    importances = sorted(
        (
            {"feature": feature, "importance": float(importance)}
            for feature, importance in zip(FEATURES, model.feature_importances_)
        ),
        key=lambda row: row["importance"],
        reverse=True,
    )
    print(format_metrics("CATH-17287 validation", validation_metrics))
    print(format_metrics("CATH-663 held-out", test_metrics))
    print("feature importance: " + ", ".join(
        f"{row['feature']}={row['importance']:.3f}" for row in importances[:10]
    ))

    if args.selection_out is not None:
        args.selection_out.parent.mkdir(parents=True, exist_ok=True)
        test_selected.to_csv(args.selection_out, index=False)
    if args.report_out is not None:
        args.report_out.parent.mkdir(parents=True, exist_ok=True)
        args.report_out.write_text(
            json.dumps(
                {
                    "features": FEATURES,
                    "validation": asdict(validation_metrics),
                    "held_out": asdict(test_metrics),
                    "feature_importance": importances,
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())