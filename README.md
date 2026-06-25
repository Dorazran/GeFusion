# GeFusion — Multi-Method Simulation Toolkit for Dopant Diffusion in SiGe Nanodevices

Kinetic Monte Carlo, global Monte Carlo sensitivity analysis (Sobol), and 2D finite-element/finite-difference modeling of how Boron, Arsenic, and Phosphorus diffuse through Silicon-Germanium (SiGe) thin films at sub-micron device scale.

This repository contains the simulation code and the raw/processed results from that work. 
This is the reproducible computational core, independent of any particular write-up.

## Background

SiGe alloys are central to modern heterojunction bipolar transistors (HBTs) and strained-channel devices because germanium incorporation tunes both the bandgap and the diffusivity of common dopants. A simple "Ge fraction changes diffusivity by some exponential factor" relationship is widely used in device design, but it is typically fit from continuum (Fick's law) models alone. This project asks a more careful question: does that phenomenological exponent actually hold up when the same physics is modeled from first principles, atom by atom?

Three independent computational methods are used to answer that:

1. **Kinetic Monte Carlo (KMC)** — a mean-field, atomistic random-walk model of individual dopant hops through a SiGe lattice, run across a 1224-condition grid (3 dopant species x 3 anneal temperatures x 17 Ge fractions x 8 independent random seeds) to extract an emergent diffusivity exponent alpha directly from simulated motion, with no continuum assumptions baked in.
2. **Global sensitivity analysis (Sobol indices)** — a variance-based Monte Carlo method that ranks which physical inputs (Arrhenius prefactor D0, activation energy Ea, temperature, Ge fraction, anneal time, and the alpha exponent itself) actually drive uncertainty in the predicted junction depth, separating what matters from what doesn't.
3. **2D finite-element/finite-difference (FEM/FD) diffusion-strain modeling** — couples the diffusion equation to mechanical strain fields in two spatial dimensions, going beyond the 1D approximations used in the continuum and KMC models.

## Key results

**KMC-emergent alpha vs. the phenomenological alpha** (from `matlab/kmc_alpha_fit_summary.csv`, 1224 underlying KMC runs):

| Species | T (C) | alpha (KMC, emergent) | alpha (phenomenological) | Difference | R-squared of fit |
|---|---|---|---|---|---|
| Arsenic | 900 | 2.57 | 2.3 | +11.7% | 0.9999 |
| Arsenic | 1000 | 2.40 | 2.3 | +4.4% | 1.0000 |
| Arsenic | 1100 | 2.23 | 2.3 | -2.9% | 1.0000 |
| Boron | 900 | -3.15 | -3.0 | -5.1% | 1.0000 |
| Boron | 1000 | -2.90 | -3.0 | +3.4% | 1.0000 |
| Boron | 1100 | -2.68 | -3.0 | +10.6% | 1.0000 |
| Phosphorus | 900 | 0.86 | 0.7 | +23.3% | 1.0000 |
| Phosphorus | 1000 | 0.80 | 0.7 | +14.6% | 1.0000 |
| Phosphorus | 1100 | 0.75 | 0.7 | +7.1% | 1.0000 |

The KMC-emergent exponents track the phenomenological values closely (within roughly 3-24%) across all nine (species, temperature) conditions, with the gap consistently narrowing at higher temperature. Phosphorus shows the largest, most systematic deviation, suggesting that its phenomenological alpha is the least well-calibrated of the three.

**Sobol global sensitivity** (first-order S_i / total-order ST_i, `matlab/sobol_results_*.csv`): for all three dopants, the Arrhenius prefactor D0 and activation energy Ea dominate the variance in predicted junction depth (D0 alone explains 55-65% of output variance), while anneal time and the alpha exponent itself contribute comparatively little on their own. This means most of the practical uncertainty in junction-depth predictions stems from how well D0 and Ea are experimentally calibrated, not from the diffusion model's structural assumptions.

Result figures:
- `matlab/sobol_indices_bar.png`, `matlab/sobol_comparison_all_dopants.png` — Sobol indices, all dopants
- `matlab/xj_distribution_histogram_*.png` — Monte Carlo output distributions per dopant
- `matlab/fem2d_concentration_field.png`, `matlab/fem2d_strain_field.png`, `matlab/fem2d_centerline_vs_1d_reference.png` — 2D FEM/FD results
- `matlab/kmc_msd_plot.png` — KMC mean-squared-displacement diagnostic
- `analysis/fig10_kmc_validation.png` — KMC-vs-phenomenological alpha comparison (the figure behind the table above)

## Repository structure

```
.
├── matlab/                       # All simulation code + results
│   ├── kmc_simulation.m          # Single-run KMC engine (reference/smoke-test)
│   ├── run_kmc_core.m            # Core KMC engine, factored for parameterized grid sweeps
│   ├── run_kmc_grid.m            # First KMC grid sweep driver (243 runs)
│   ├── run_kmc_grid_v2.m         # Deeper KMC grid sweep driver (1224 runs, used for results above)
│   ├── global_mc_sobol.m         # Global Monte Carlo + Sobol sensitivity analysis
│   ├── fem_2d_diffusion_strain.m # 2D FEM/FD diffusion-strain solver
│   ├── replot_sobol_bar.m        # Re-plotting utility for Sobol bar charts
│   ├── kmc_grid_v2_runs/         # Raw per-run checkpoints (1224 runs, full provenance)
│   ├── kmc_grid_v2_aggregated.csv, kmc_grid_v2_results.csv, kmc_alpha_fit_summary.csv
│   ├── sobol_results_*.csv/.mat, fem2d_results.mat
│   └── *.png                     # Result figures
├── analysis/
│   ├── regen_fig10.py            # Builds the KMC-vs-phenomenological comparison figure from the CSVs above
│   └── fig10_kmc_validation.png
├── LICENSE
└── README.md
```

## Code excerpt

The KMC core engine (`matlab/run_kmc_core.m`) takes a parameter struct so the same engine can be called across the full grid with independent RNG seeds per replicate:

```matlab
function results = run_kmc_core(opts)
% Core mean-field KMC engine. INPUT: opts struct with fields:
%   species              'Boron' | 'Arsenic' | 'Phosphorus'
%   X_GE                 uniform Ge fraction, 0<=X_GE<=1
%   T_ANNEAL_K           anneal temperature, K
%   rng_seed             integer seed -- REQUIRED for genuine statistical
%                        independence between replicate runs
%   N_DOPANTS, NX, NY, NZ, MIN_HOPS_PER_DOPANT, CONVERGENCE_TOL, ...
%
% OUTPUT: results struct with D_eff_final, D_SiGe_analytic_cm2_s,
%   relative_difference_pct, MSD_vs_t_*, params, elapsed_wallclock_sec
```

## Running the simulations

Requires MATLAB (no special toolboxes beyond base MATLAB for the KMC/Sobol/FEM scripts). From the `matlab/` directory:

```matlab
kmc_simulation          % single smoke-test KMC run
run_kmc_grid_v2         % full 1224-run grid (multi-hour wall-clock)
global_mc_sobol         % Sobol global sensitivity analysis
fem_2d_diffusion_strain % 2D FEM/FD diffusion-strain model
```

`analysis/regen_fig10.py` (Python, requires numpy/matplotlib) rebuilds the comparison figure directly from the result CSVs in `matlab/`.

## Credit

This work was carried out as part of a microelectronics course project under the guidance of **Prof. Gadi Golan**, Department of Electrical and Electronics Engineering, **Ariel University**.

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Dor Azran, Ariel University.
