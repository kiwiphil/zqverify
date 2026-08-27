# Theorem-to-code map

| Manuscript object | Code |
|---|---|
| finite active character group / quotient | `groups.py` |
| signed quotient word length | `groups.signed_word_length` |
| detector coefficient and optimization | `detector.py` |
| cubical incidence conventions | `lattice.py` |
| signed multicolour current and source | `currents.py` |
| deterministic quotient extraction | `currents.extract` |
| positive-character finite Fourier transform | `fourier.py` |
| tiny-volume flux sum | `flux.py` |
| detector catalogue | `search.py`, `examples/search_small_groups.py` |
| direct Wilson-vs-detector audit | `examples/finite_volume_bound_audit.py` |

The analytic source-sector domination theorem is checked only through finite Wilson ratios; the program does not prove it for arbitrary volumes. The projected-filling theorem and history-tree estimate likewise remain analytic theorems.
