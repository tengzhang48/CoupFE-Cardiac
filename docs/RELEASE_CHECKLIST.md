# Release checks

Run these checks from a clean clone of the proposed public revision.

## Dependency and repository

- [ ] `pyproject.toml` pins CoupFE to one full commit on the public HTTPS
  repository.
- [ ] The pinned Core commit is anonymously reachable and the installed
  `direct_url.json` passes `.github/scripts/check_runtime_core.py`.
- [ ] The release commit has reviewed history and author metadata, and the
  release worktree is clean.
- [ ] `git diff --check` and `git fsck --full --strict` pass.

## License and credit

- [ ] `LICENSE`, `NOTICE`, `THIRD_PARTY_NOTICES.md`, `docs/LICENSE.md`, and both
  files under `LICENSES/` are present in source and built artifacts.
- [ ] `activation.py`, `fiber_crosscheck.py`, and `structural_directions.py`
      retain their SPDX identifiers,
  immutable upstream revisions, file-level boundaries, and author credit.
- [ ] The Finsberg/Sundnes/van den Brink activation port uses the exact Zenodo
  v1.0.0 deposit's CC BY 4.0 record, names its archive and hashes, links the
  license, and indicates the changes made by CoupFE-Cardiac. No downstream-
  created MIT copyright notice is presented as an upstream license.
- [ ] The benchmark paper and Zenodo dataset remain cited; the raw archive and
  team pickles stay unbundled, while retained team statistics, envelopes,
  selected curves, figures, and RED values name the creators and DOI, link
  CC BY 4.0, and explicitly indicate CoupFE-Cardiac's transformations.

## Executable gates

- [ ] `python -m pytest -q` passes with no optional-gate skip counted as a
  success.
- [ ] `python -m pytest -q -m slow` passes with `gfortran` available.
- [ ] The slow gate completes all 500 current Case A steps and checks the
  Hex8-sampled retained time history and peak field within the recorded
  tolerances.
- [ ] `python -m pytest -q -m mpi` passes in the documented matched PETSc/MPI
  environment; retain the 1/2/4-rank output for every script.
- [ ] Before closed-mesh MPI refinement, retain serial and 1/2/4-rank results
  across the 0.20--0.32 s snap window for the same pointwise-`kappa`,
  consistent-mass, Laplace/GP-direct configuration. Each shortened run records
  `load_horizon=1.0 s`, and its load arrays equal the production prefix. Run
  `compare_mpi_rank_gate.py`; do not substitute manual log inspection or a
  smaller smoke mesh for this gate.
- [ ] The Q1/P0 invalid-trial controls pass: Core exercises the operator
  `max_step` bound; a real PETSc `bt` check recovers only
  `InvalidDeformationError` from the residual callback and records the
  rejection; unrelated residual and Jacobian exceptions remain fatal.
- [ ] Any closed-domain Case B claim passes the pre-solve identity gate before
  time integration: units/extents, every exterior face classified exactly
  once, positive Gauss and extended reference Jacobians, boundary measures,
  pressure orientation/resultant/moment, and Robin symmetry. Retain the audit
  JSON with the run.
- [ ] A paper-parameter reproduction records `kappa=1.0e6 Pa` and
  `eta=100 Pa s`; stopped or failed changed-parameter sensitivity controls are
  not promoted to benchmark results.
- [ ] No Case A or Case B numerical claim is made without a retained result,
  environment, command, external-data identity, and RED output from its exact
  clean application checkpoint and approved Core revision.
- [ ] A retained direct FEniCS comparison verifies expected hashes for all five
  input roles, accepts only the corrected material-law ID and fixed closed Case
  B setup, records serial or validated MPI provenance, and uses the fixed
  common grid, snap window, and onset definition.

## Artifacts and public boundary

- [ ] `python -m build --sdist --wheel` succeeds from the clean revision.
- [ ] `python -m twine check --strict dist/*` passes.
- [ ] `python .github/scripts/check_release_artifacts.py dist --source-root .`
  passes without an audit override and reports the exact reviewed inventories.
- [ ] The source archive contains the runnable application, docs, and tests.
  The wheel is intentionally metadata/notices-only and is not described as an
  application-code distribution.
- [ ] Every retained report JSON and normalized transcript is present in the
  source tree and sdist; its provenance, source hashes, semantics, and bounded
  interpretation pass the release guard.
- [ ] The supplied FEniCS point-stress arrays are not used as a quantitative
  oracle. A stress comparison, if included, reloads accepted `u/v/a`, uses a
  corrected stress invariant, declares element-interior/quadrature locations,
  and records projection and physical sample distance.
- [ ] A dedicated scan finds no credentials, private coordination material,
  machine-local paths, generated kernels, shared libraries, generated binary
  result archives, plots, raw external data, or unreviewed files in the public
  commit and artifacts.

## Evidence record

- [ ] Retain UTC environment/package/toolchain versions, exact commands and
  exit codes, app/Core SHAs and tree states, source/artifact inventories,
  artifact hashes, and the final remote ref inspection.
- [ ] Record optional gates that were not run as not run; do not convert them to
  passes.
- [ ] Keep the README Code, Tests, Examples, Retained results, and Limitations
  sections consistent with the retained evidence and the geometry,
  fiber-sampling, time-integration, Case B, and MPI boundaries.

## Open blockers

Items known to be unresolved. Each must be closed or consciously accepted and
recorded before the public release is tagged.

### B1. Closed benchmark geometry is the CLI default — resolved

**Status: closed 2026-08-12.** Both drivers now default to
`closed-multiblock`, and selecting `polar-ring` with a positive apex
offset prints an explicit not-the-benchmark-geometry warning. The
polar-ring option is retained for archived-evidence reproduction; every test
that exercises it passes `--mesh-topology polar-ring` explicitly. This preserves
the historical records without exposing the non-benchmark open-tip geometry as
an implicit user choice.

**Related, lower severity.** `post.py` silently infers `polar_ring` for any
archive lacking a `mesh_topology` field. That behaviour is deliberate and
load-bearing for reading archived NPZs, and is guarded by a test, but it should
record that the topology was *inferred* rather than read.
