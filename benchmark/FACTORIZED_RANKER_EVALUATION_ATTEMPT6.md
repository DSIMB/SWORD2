# Factorized structural ranker: locked attempt 6

Date: 2026-08-13

## Outcome

Attempt 6 is a complete, valid measured failure. The frozen factorized ranker
passed 2 of 9 predeclared gates and failed 7. It must remain experimental and
opt-in; plain SWORD2 invocations continue to use the legacy selector.

The factorized model made actual learned selections for 642 of 663 CATH chains
(96.83257918552036% structural coverage). The other 21 chains were predeclared
input-only abstentions because at least one candidate residue lacked complete,
finite N/CA/C/O backbone evidence. Their legacy answers remain available for
usability but are not counted as factorized predictions.

No benchmark command or metric evaluation was rerun after the valid evaluator
exit. The earlier invalid attempts remain preserved as audit evidence.

## Predeclared gate results

All performance values below are conditional on the exact 642-chain
structurally eligible set.

| Gate | Exact value | Rule | n | Result |
|---|---:|:---:|---:|:---:|
| overall NDO | 0.7980101951659153 | > 0.8389 | 642 | FAIL |
| factorized-minus-Merizo NDO CI low | -0.059402207277476755 | > 0.0 | 642 | FAIL |
| factorized-minus-Chainsaw NDO CI low | 0.058115974018791136 | > 0.0 | 602 | PASS |
| domain-count accuracy | 0.7024922118380063 | >= 0.745 | 642 | FAIL |
| boundary F1 at 10 residues | 0.5894570815131562 | >= 0.62 | 642 | FAIL |
| contiguous factorized-minus-standalone NDO | -0.01679479417379004 | >= -0.005 | 381 | FAIL |
| discontinuous factorized-minus-standalone NDO | -0.01580869266609619 | >= -0.005 | 261 | FAIL |
| median paired runtime ratio | 1.2851294463681646 | <= 1.15 | 642 | FAIL |
| ratio of peak RSS maxima | 0.9999509539457551 | <= 1.1 | 642 | PASS |

The 10,000-draw seed-37 paired bootstrap estimated:

- factorized minus Merizo mean NDO: -0.04277394245210569, 95% CI
  [-0.059402207277476755, -0.026460949829754114];
- factorized minus Chainsaw mean NDO: 0.07813088211148893, 95% CI
  [0.058115974018791136, 0.0986477460699529].

Conditional mean NDO was 0.7980101951659153 for factorized,
0.7780834731887176 for runtime legacy, 0.8407841376180212 for Merizo, and
0.7249673276681073 for Chainsaw on its 602-chain eligible intersection.
Thus the model improves over its runtime legacy selector by about 0.01993 NDO
on eligible chains, but is materially below Merizo and does not meet the locked
quality, boundary, count, cohort, or runtime criteria.

Resource diagnostics were a runtime-ratio p95 of 1.5534250035390835,
factorized maximum RSS of 81,552 kB, and legacy maximum RSS of 81,556 kB.

## Coverage and denominator integrity

- Full audit denominator: 663 chains.
- Factorized eligible and scored: 642.
- Predeclared structural abstentions: 21, all
  `incomplete_backbone`.
- Factorized candidate exclusions after ranking began: 0.
- Eligible paired resource measurements: 642.
- Merizo full collection: 663.
- Frozen Chainsaw successes: 623; expected failures: 40.
- Chainsaw comparison after structural eligibility intersection: 602.
- Resource order SHA-256:
  `4ce265d73f087bc5a045492921a3fe0a1686992cfe9aaa9a2a69e3f46743f933`.

This policy does not declare an incomplete biological structure unusable.
It only abstains from applying frozen model v1 to a representation it was not
trained or calibrated to interpret. A gap-aware model would require separate
training, calibration, missingness features, and untouched evaluation data.

## Rollout decision

The selector remains opt-in for three independent reasons:

1. seven of nine performance/resource gates failed;
2. learned structural coverage is below 100%;
3. the frozen cache contract is not default-promotion compatible
   (`whole_chain_legacy_fallback`, no typed-context reload).

The exact invocation remains:

```bash
sword2 -i structure.pdb -o results --use-factorized-ranker
```

Omitting `--use-factorized-ranker` remains the plain legacy escape hatch.

CATH-663 informed earlier diagnosis and protocol design. These results are
therefore engineering evidence, not a publication-quality generalization
claim; an untouched external benchmark is required for that.

## Audit authority

- Model manifest SHA-256:
  `b45f6d311a33a790cfe03615f347a4626e5c0a6412309c39e4d07882f98a09c0`.
- Runtime manifest SHA-256:
  `770c3d348365b7f744f2eaee964b61cc57615c10f8ca9a47a061ac03ee1c343e`.
- Locked binary SHA-256:
  `0d9e5fdba57c30b9b0e238876011815027bab365040d0c605f74fba1987f67c5`.
- Attempt-6 intent SHA-256:
  `b387e92bf551bb12b60398123dc98b0cc76c3e82188439d561d67a48a2d5971e`.
- Eligibility manifest SHA-256:
  `658d37a987eb0aeac9e13f6c6a1574956c18b0a125cfc32fbc7a84d668eb7558`.
- Coverage attestation SHA-256:
  `9956de5a07984229a20ea02029918b41d4fa71cc653fa3a13751ad05dae0bd9b`.
- Acceptance JSON SHA-256:
  `654c0bc96223f5ab21391db735d71e7da8c8cf958011a55477f6533d12b3e282`.
- Acceptance Markdown SHA-256:
  `530c46ed3c1a30c6721f953dfde49ca43756ad6a687a847eb3617c950793ef49`.
- Frozen event ledger SHA-256:
  `167f10956fe1cb1f9a4e66e7d42d35c17e848f96552afa5d5bf2713d84219976`.
- Execution receipt SHA-256:
  `c357a09d7e75ac7296736c8b91efef5bcf527c00ad3b001556c05653abe07dee`.

The receipt is canonical and independently checked against the immutable
16-event chain, five normalized commands and logs, exit codes, three result
manifests, eleven evidence files, coverage, and acceptance artifacts.

## Formal promotion-validator limitation

The frozen validator's commit-binding layer hardcodes the old generic receipt
and report paths. Attempt 6 deliberately uses collision-safe, attempt-specific
paths after preserving attempts 3 through 5, while
`benchmark/REPORT.md` was already an unrelated untracked user file and was
left byte-for-byte untouched. Consequently, no valid `PROMOTE` token can be
issued from this checkout without either mutating user-owned data or changing
the verifier after observing results. Both would be scientifically unsafe.

This packaging limitation cannot turn failed gates into a promotion. The
fail-closed operational outcome remains to retain opt-in. A future administrative
repair may add attempt-path parameters to a newly frozen validator and rebind
these exact existing bytes; it must not rerun or rescore CATH-663.
