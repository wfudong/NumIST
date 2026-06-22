# Project Inventory

## Existing Root-Level Files

- `README.md`: short description for the paper-reproduction workflow.
- `paper_examples_kdv.jl`: headless script for Section 9.3 of `152.pdf`.
  It mixes reusable ETDRK4 numerics, paper initial data, plotting, and artifact
  writing in one file.
- `PaperExamples_KdV.ipynb`: notebook wrapper around `paper_examples_kdv.jl`,
  plus custom evolution/animation cells.
- `ISTM_finite_section_KdV.ipynb`: standalone finite-section inverse-scattering
  reconstruction code. This is not yet migrated into `NumISTNext`.
- `KdVSolitons.ipynb`: exploratory Crum/Darboux examples using
  `KdVCrumProject`.
- `CrumSolitonGasInteraction.ipynb`: combines Crum-generated soliton-gas
  initial data with the ETDRK4 KdV solver.
- `figures/`, `soliton_gas_evolution.gif`, and notebook-local GIFs: generated
  output currently tracked beside source files.

## Existing `KdVCrumProject`

- `Project.toml` and `Manifest.toml`: package-like Julia environment.
- `src/CrumTypes.jl`: problem/result types and spectrum validation.
- `src/SeedTrajectories.jl`: sampled seed trajectories and stable rescaling.
- `src/Crum1Fold.jl`, `src/Crum2Fold.jl`: stable branch formulas and dressing
  maps.
- `src/CrumAlgorithm.jl`: paired Crum transform and naive Darboux reference.
- `src/Utilities.jl`: exact one-soliton, norming-constant evolution, invariant
  errors, spectra helpers, and clipping.
- `src/Visualization.jl`: plotting helpers; this makes the old main module load
  plotting packages even for pure numerical use.
- `notebooks/`: one-, two-, many-soliton, and interactive demo notebooks.

## Rebuild Decisions

- Keep all legacy files in place until the new package has matching tests and
  examples.
- Make `NumISTNext/src` the new source of truth for reusable code.
- Keep generated files under `NumISTNext/output/` or top-level `figures/`, never
  mixed into source directories.
- Treat notebooks as consumers of the package, not as primary implementation
  storage.
- Migrate finite-section IST after the KdV solver and Crum workflows are stable.
