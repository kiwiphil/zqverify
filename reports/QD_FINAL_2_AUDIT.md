# QD-FINAL.2 release audit

## Verdict

**PASS as a local public-release candidate, pending creation of the actual persistent public repository/archive identifier.**

## Manuscript synchronization

The v6 manuscript and this release use the same conventions:

- one representative of every non-self-inverse active-character pair;
- arbitrary subgroup detector kernels;
- visible colours are those not in the kernel;
- extraction requires nonzero projected incidence at the first defective edge but does not require monotone defect decrease;
- source-sector domination and the history estimate remain analytic proofs, not computational assumptions.

The release reconstruction uncovered and corrected one convention inconsistency in the v5 manuscript: its Z6 example simultaneously listed colours 2 and 4=-2. The v6 example now consolidates that inverse pair exactly as required by the setup.

## Tests

- Source-tree test suite: 76/76 passing.
- Built wheel: successful.
- Wheel installed into an isolated target: successful.
- Same 76 tests run against the installed wheel: passing.

## Catalogue audit

Fresh catalogue: 9,156 detector problems across Z2, Z3, Z4, Z5, Z6, Z2^2, Z2^3, Z2 x Z4, and Z3^2, using every generating canonical active-colour subset of size at most 3 and effective activities in {0.03, 0.12, 0.55}.

Outcomes:

- strict noncyclic-detector advantage: 1,960;
- crude total rho >= 1 but an admissible quotient detector exists: 2,109.

## Direct finite-volume audit

The historical 214-case family is reproduced exactly:

- 214/214 inequalities pass;
- 162 product-group one-plaquette d=2 cases;
- 131 have a nontrivial optimizing kernel;
- 31 have a noncyclic optimizing quotient;
- maximum exact Wilson expectation / detector upper bound = 0.4849636282734537.

## Proof boundary

These computations do not prove source-sector domination, the general projected rectangular filling theorem, the infinite history-tree bound, or thermodynamic-limit inheritance. Those results remain analytic theorems in the manuscript.

## Remaining release action

Before arXiv/journal submission, create the actual public repository or archival deposit and replace the manuscript's placeholder sentence with its persistent URL/DOI. No public identifier is fabricated here.
