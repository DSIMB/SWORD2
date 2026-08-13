# Factorized structural ranker acceptance

Structural coverage: 642/663 (96.832579%); abstentions: 21.

All performance gates below are conditional on structurally eligible chains. Legacy fallbacks for abstained chains are not factorized scores.

| Gate | Value | Rule | n | Result |
|---|---:|:---:|---:|:---:|
| overall_ndo | 0.798010195166 | > 0.8389 | 642 | FAIL |
| merizo_ndo_ci_low | -0.0594022072775 | > 0 | 642 | FAIL |
| chainsaw_ndo_ci_low | 0.0581159740188 | > 0 | 602 | PASS |
| domain_count_accuracy | 0.702492211838 | >= 0.745 | 642 | FAIL |
| boundary_f1_10 | 0.589457081513 | >= 0.62 | 642 | FAIL |
| contiguous_standalone_ndo_delta | -0.0167947941738 | >= -0.005 | 381 | FAIL |
| discontinuous_standalone_ndo_delta | -0.0158086926661 | >= -0.005 | 261 | FAIL |
| runtime_median_ratio | 1.28512944637 | <= 1.15 | 642 | FAIL |
| rss_maxima_ratio | 0.999950953946 | <= 1.1 | 642 | PASS |

All gates measured: true
All gates pass: false
