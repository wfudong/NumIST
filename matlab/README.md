# MATLAB realization of `152.pdf`, Section 9

Run from MATLAB:

```matlab
cd('C:\Users\WFD\Desktop\Research Code\numerical IST\NumIST\matlab')
run_section9_examples
```

For experimentation, use:

```matlab
main_kdv_playground
```

Edit the first block of `main_kdv_playground.m` to choose the initial data,
the target time or times, the plotting window, and the time step.

For animations of the three paper examples, use:

```matlab
animate_section9_examples
```

This opens live MATLAB animations and saves GIFs to `animations_section9/`.
To animate just one example:

```matlab
animate_section9_examples('Example', 'example93')
```

The script writes PNG reproductions of the Section 9 examples to
`figures_section9/`, and also writes PDF copies when the local MATLAB
graphics/font setup supports vector export:

- `figure1_example91`: `q(x) = x exp(-x^2)`, `t = 0, 0.08, 0.16`
- `figure2_example92`: exact one-soliton with `c = pi`, `t = 0, 0.5, 1`
- `figure3_example93_initial`: the piecewise initial profile
- `figure4_example93`: Example 9.3 evolved to `t = 0, 0.015, 0.03`

The helper `jost_i_over_2_spps.m` implements the SPPS initial-value step
described in Section 9.1 of `152.pdf`, using the construction from
`sppsmethod.pdf`.

To reproduce the figures with the actual ISTM finite-section method from
`152.pdf` instead of ETDRK time stepping, run:

```matlab
run_section9_istm_paper_method
```

This writes figures and diagnostics to `figures_section9_istm/`.

To put the ISTM and ETDRK results on the same plots and also plot their
pointwise differences, run:

```matlab
compare_section9_methods
```

This writes overlays, difference plots, and a comparison summary to
`figures_method_comparison/`.

To compare the runtime of the paper ISTM method against ETDRK4 without
including figure export time, run:

```matlab
benchmark_section9_methods
```

This writes a timing table and log-scale bar plot to
`figures_speed_comparison/`.

For a faster paper-style figure-resolution timing, for example with 501
x-points on each displayed interval, use:

```matlab
benchmark_section9_methods('NumTheta', 10000, 'NumX', 501, ...
    'OutputDir', 'figures_speed_comparison_paper_style')
```

To animate the ISTM evolution of the reflectionless initial data
`q(x) = -110 sech^2(x)`, run:

```matlab
animate_sech110_istm
```

This uses the exact ten-soliton scattering data and writes a smooth MP4,
waterfall plot, MAT data file, and summary to `animations_sech110_istm/`.
Add `'SaveGif', true` if a GIF copy is needed.

For an extra-smooth long render of the same large window, use:

```matlab
animate_sech110_istm('SaveMat', false, 'MakeWaterfall', false, ...
    'NumX', 16000, 'NumFrames', 1200, 'FrameRate', 60, ...
    'VideoQuality', 100, ...
    'OutputDir', 'animations_sech110_istm_ultra_smooth')
```
