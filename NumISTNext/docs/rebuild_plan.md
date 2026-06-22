# Rebuild Plan

## Phase 1: Safe Package Skeleton

- Add `NumISTNext/` beside the legacy code.
- Separate reusable KdV spectral evolution from figure-generation scripts.
- Implement Section 9 finite-section inverse scattering and Section 4
  scattering-data time evolution as the primary KdV evolution path.
- Move the Crum numerical core without importing plotting packages from the main
  module.
- Keep plotting dependencies in `scripts/Project.toml`, separate from the core
  package environment.
- Add `notebooks/playground.ipynb` as the first clean interactive notebook.
- Add small tests for the one-soliton Crum case and KdV zero solution.

## Phase 2: Reproducible Examples

- Convert `paper_examples_kdv.jl` into a script that calls `NumISTNext`.
- Move custom animation and soliton-gas workflows from notebooks into scripts.
- Keep notebooks short: activate the project, call package/script functions, and
  display results.

## Phase 3: Finite-Section IST

- Extract the Section 9.1 direct-scattering workflow for computing scattering
  data from an initial potential. The inverse finite-section step and
  scattering-data time evolution already live in `src/ist/finite_section.jl`.
- Add tests against Example 9.1 and Example 9.3 once direct reflection
  coefficients are produced by the package rather than supplied externally.

## Phase 4: Cleanup Gate

- Compare generated figures and validation summaries against the legacy outputs.
- Add a top-level migration note once the new commands reproduce the old public
  artifacts.
- Only then decide whether to archive or remove obsolete legacy notebooks.
