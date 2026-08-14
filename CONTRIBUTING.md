# Contributing to CoupFE-Cardiac

This repository is an application of the separately maintained CoupFE core. Do not copy or edit
core modules here. The public dependency is pinned in `pyproject.toml`; a local editable core may
be substituted deliberately during development.

Keep the application boundary narrow:

- Cardiac geometry, constitutive configuration, activation, boundary operators, diagnostics, and
  drivers live under `examples/`.
- Cardiac owns ventricular mesh generation, facet/region interpretation, and
  fiber/sheet attachment. Serial operators receive the resulting neutral NumPy
  coordinate/connectivity arrays; distributed drivers wrap those arrays in
  Core's `KernelMeshView` for partitioning. Do not move these adapters or
  benchmark semantics into Core.
- Generic assembly, code generation, compiled-element runtime, and PETSc/MPI infrastructure belong
  in CoupFE core and need their own tests there.
- Run the fast cardiac tests before every change. Compiled-kernel and MPI gates are explicit
  opt-ins; never interpret a skipped optional gate as a pass.
- Never commit generated kernels, shared libraries, result archives, plots, reference data,
  machine-local paths, credentials, or private coordination notes.

Keep evidence claims within the executable record:

- The public default is backward Euler. Newmark kinematics are an experimental
  option until the material-rate and boundary-rate discretizations are assessed
  together.
- The procedural demo uses a small open-apex mesh to avoid collapsed elements.
  It is not a direct reproduction of the closed benchmark geometry.
- Four MPI drivers compare distributed and serial implementations of the same
  small problems within a stated tolerance; the scaling driver compares with
  its distributed rank-1 result. They do not independently validate cardiac
  mechanics or establish general performance.
- External benchmark comparisons require a user-supplied dataset and a result
  file whose application and Core provenance pass the comparison tool's checks.

When reporting a result, retain the exact command, app and core revisions, environment, result
checksum, and relevant console output. See `docs/BENCHMARK_COMPARISON.md` and
`examples/REFERENCES.md`.
