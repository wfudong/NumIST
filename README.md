# NumIST

Julia notebooks and scripts for the numerical examples in `152.pdf`, the paper
"Realization of the inverse scattering transform method for the Korteweg-de Vries
equation."

## Files

- `ISTM_finite_section_KdV.ipynb`: finite-section inverse scattering
  reconstruction for the KdV equation, with a reflectionless soliton validation.
- `paper_examples_kdv.jl`: headless Julia script that regenerates the figures
  for all examples in Section 9.3 of the paper.
- `PaperExamples_KdV.ipynb`: notebook version of the full example reproduction.
- `figures/`: saved PNG and PDF figures matching Figures 1-4 in the paper,
  plus `validation_summary.txt`.

## Regenerate the figures

From the repository root:

```bash
julia --startup-file=no paper_examples_kdv.jl
```

This writes:

- `figures/figure1_example91.{png,pdf}`
- `figures/figure2_example92.{png,pdf}`
- `figures/figure3_example93_initial.{png,pdf}`
- `figures/figure4_example93.{png,pdf}`
- `figures/validation_summary.txt`

The script uses `FFTW` and `Plots`. The notebook can also be executed with
Jupyter using the Julia kernel.
