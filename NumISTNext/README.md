# NumISTNext

This folder is the additive rebuild of the original `NumIST` workspace. The old
code, notebooks, generated figures, and PDFs stay untouched at the repository
root and in `KdVCrumProject/`.

## Layout

- `src/kdv/`: reusable KdV time-evolution and paper-example initial data.
- `src/crum/`: migrated Crum/Darboux reflectionless-potential core, without the
  old plotting dependency in the main module.
- `scripts/`: reproducible commands that write generated artifacts under
  `output/`.
- `notebooks/`: interactive Julia notebooks for exploration.
- `test/`: small behavioral checks to protect the migration.
- `docs/`: project inventory and rebuild plan.

## Commands

From the repository root:

```bash
julia --project=NumISTNext -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=NumISTNext/scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=NumISTNext/scripts NumISTNext/scripts/reproduce_example92_ist.jl
```

The last command writes the Section 9 finite-section IST reproduction of paper
Example 9.2 to `NumISTNext/output/figures_ist/`.

Open `notebooks/playground.ipynb` when you want an interactive place to change
parameters and plot soliton or KdV evolution examples.
