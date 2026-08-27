# zqverify 1.0.0

Executable verification companion for **Quotient detectors and strong-coupling area laws in finite Abelian lattice gauge theory**.

The package mirrors finite objects used in the manuscript: products of cyclic groups, subgroup quotients, signed quotient word metrics, detector optimization, oriented cubical incidence, multicolour signed currents, the deterministic quotient-source extraction, finite Fourier coefficients, and exact tiny-volume Wilson expectations.

## Scope

The software is a regression, convention, and finite-instance verification layer. It is **not** a proof assistant and does not replace the analytic proofs of source-sector domination, the general rectangular filling theorem, the history-tree bound, or the thermodynamic-limit statements.

## Reproduce

```bash
python -m pip install -e '.[test]'
pytest -q
python examples/search_small_groups.py
python examples/finite_volume_bound_audit.py
```

The finite-volume audit deliberately reproduces the 214-case family described in the manuscript: 162 product-group one-plaquette cases in d=2, 48 cyclic d=2 rectangle cases, and 4 Z2 one-cube d=3 cases.

See `docs/THEOREM_TO_CODE.md` and `docs/COMPUTATIONAL_SCOPE.md`.
