import numpy as np
import pandas as pd

from benchmark.evaluate_structural_ranker import metrics, select_top1


def test_select_top1_picks_highest_score_per_chain():
    frame = pd.DataFrame(
        {
            "chain_id": ["a", "a", "b", "b"],
            "ndo": [0.2, 0.9, 0.7, 0.4],
            "d_count_acc": [0.0, 1.0, 1.0, 0.0],
            "boundary_f1_10": [0.0, 1.0, 1.0, 0.0],
            "matched_dice": [0.2, 0.9, 0.8, 0.4],
            "n_pred_domains": [1, 2, 2, 3],
        }
    )

    selected = select_top1(frame, np.array([0.1, 0.9, 0.8, 0.2]))

    assert selected["ndo"].tolist() == [0.9, 0.7]
    assert metrics(selected).n_chains == 2
    assert metrics(selected).ndo == 0.8