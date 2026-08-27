# Reproducibility

Tested release: `zqverify 1.0.0`.

```bash
python -m pip install -e '.[test]'
python -m pytest -q
python examples/search_small_groups.py
python examples/finite_volume_bound_audit.py
```

Release results:

- 76/76 unit and parametrized regression tests pass.
- Small-group detector catalogue: 9,156 problems; 1,960 strict noncyclic-detector advantages; 2,109 cases with crude total activity at least one but an admissible quotient detector.
- Direct finite-volume audit: 214/214 detector inequalities pass; 162 are quotient-sensitive product-group one-plaquette cases; 131 have a nontrivial optimizing kernel and 31 a noncyclic optimizing quotient; maximum observed exact-Wilson/bound ratio is 0.4849636282734537.
- The built wheel was installed into an isolated target and the same 76 tests passed against the wheel.

All computations are finite-instance checks. The analytic manuscript remains the proof of the general results.
