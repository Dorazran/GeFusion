# -----------------------------------------------------------------------------
# Copyright (c) 2026 Dor Azran, Ariel University
# Licensed under the MIT License. See LICENSE file in the project root.
# -----------------------------------------------------------------------------

"""
Regenerate Fig. 4 (fig10_kmc_validation.png) for the AESMT manuscript.
Left panel : ln(D_eff) vs x_Ge per species at 1000 C (T=1273.15 K), with linear KMC fits.
Right panel: species-averaged alpha_KMC vs alpha_phenom (bar comparison).
Data sources:
  - kmc_grid_v2_aggregated.csv   (D_eff_mean_cm2_s per species/T/x_Ge, for left panel + linear fit)
  - kmc_alpha_fit_summary.csv    (alpha_KMC / alpha_phenom per species/T, for right panel)
Rebuilt with much larger fonts/figure size so it is legible inside the AESMT
template's narrow two-column body (~3.0 in wide).
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MATLAB_DIR = "/Users/dorazran/Desktop/diffusion_SiGe_paper/revision_v3_deep/matlab"
PAPER_DIR = "/Users/dorazran/Desktop/diffusion_SiGe_paper/revision_v3_deep/paper"

GRID_CSV = os.path.join(MATLAB_DIR, "kmc_grid_v2_aggregated.csv")
ALPHA_CSV = os.path.join(MATLAB_DIR, "kmc_alpha_fit_summary.csv")
OUT_PNG = os.path.join(PAPER_DIR, "fig10_kmc_validation.png")

T_TARGET_K = 1273.15  # 1000 C

# -------------------- load grid data (left panel) --------------------
species_data = {}  # species -> (x_Ge array, lnD array)
with open(GRID_CSV, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        T_K = float(row["T_K"])
        if abs(T_K - T_TARGET_K) > 0.01:
            continue
        sp = row["species"]
        x_Ge = float(row["x_Ge"])
        D = float(row["D_eff_mean_cm2_s"])
        species_data.setdefault(sp, []).append((x_Ge, math.log(D)))

for sp in species_data:
    species_data[sp].sort(key=lambda t: t[0])

# -------------------- load alpha summary (right panel) --------------------
alpha_rows = []
with open(ALPHA_CSV, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        alpha_rows.append(row)

species_order = ["Boron", "Arsenic", "Phosphorus"]
avg_alpha_kmc = {}
avg_alpha_phenom = {}
for sp in species_order:
    vals_kmc = [float(row["alpha_KMC"]) for row in alpha_rows if row["species"] == sp]
    vals_phenom = [float(row["alpha_phenom"]) for row in alpha_rows if row["species"] == sp]
    avg_alpha_kmc[sp] = sum(vals_kmc) / len(vals_kmc)
    avg_alpha_phenom[sp] = sum(vals_phenom) / len(vals_phenom)

# -------------------- plotting --------------------
# Stacked (one panel below the other) layout, much larger fonts, so the
# figure stays legible after being shrunk to the AESMT two-column body width.
plt.rcParams.update({
    "font.size": 22,
    "axes.titlesize": 24,
    "axes.labelsize": 24,
    "xtick.labelsize": 21,
    "ytick.labelsize": 21,
    "legend.fontsize": 19,
    "axes.linewidth": 1.6,
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 17))

colors = {"Boron": "#1f77b4", "Arsenic": "#d62728", "Phosphorus": "#2ca02c"}
markers = {"Boron": "o", "Arsenic": "s", "Phosphorus": "^"}

for sp in species_order:
    pts = species_data.get(sp, [])
    if not pts:
        continue
    xs = np.array([p[0] for p in pts])
    ys = np.array([p[1] for p in pts])
    ax1.plot(xs, ys, markers[sp], color=colors[sp], markersize=13,
              markeredgecolor="black", markeredgewidth=1.0, label=f"{sp} (KMC)")
    # linear fit: ln(D) = ln(D0) + alpha * x_Ge
    slope, intercept = np.polyfit(xs, ys, 1)
    xfit = np.linspace(xs.min(), xs.max(), 50)
    ax1.plot(xfit, slope * xfit + intercept, "-", color=colors[sp], linewidth=3.2,
              label=f"{sp} fit (alpha={slope:.2f})")

ax1.set_xlabel(r"Ge fraction $x_{Ge}$", fontweight="bold")
ax1.set_ylabel(r"$\ln(D_{eff})$  [cm$^2$/s]", fontweight="bold")
ax1.set_title("KMC-derived ln($D_{eff}$) vs. $x_{Ge}$ at 1000 C", fontsize=24)
ax1.legend(loc="best", framealpha=0.95, edgecolor="black")
ax1.grid(True, linewidth=0.8, alpha=0.6)
ax1.tick_params(width=1.5, length=8)

# -------------------- bottom panel: bar comparison --------------------
x = np.arange(len(species_order))
width = 0.35
vals_kmc = [avg_alpha_kmc[sp] for sp in species_order]
vals_phenom = [avg_alpha_phenom[sp] for sp in species_order]

ax2.bar(x - width / 2, vals_kmc, width, label=r"$\alpha_{KMC}$ (avg.)",
        color="#4c72b0", edgecolor="black", linewidth=1.2)
ax2.bar(x + width / 2, vals_phenom, width, label=r"$\alpha_{phenom}$",
        color="#dd8452", edgecolor="black", linewidth=1.2)
ax2.axhline(0, color="black", linewidth=1.2)
ax2.set_xticks(x)
ax2.set_xticklabels(species_order, fontsize=21)
ax2.set_ylabel(r"Ge-enhancement exponent $\alpha$", fontweight="bold")
ax2.set_title("Species-averaged alpha_KMC vs. alpha_phenom", fontsize=24)
ax2.legend(loc="best", framealpha=0.95, edgecolor="black")
ax2.grid(True, axis="y", linewidth=0.8, alpha=0.6)
ax2.tick_params(width=1.5, length=8)

fig.tight_layout(pad=2.5)
fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight")
print("Saved:", OUT_PNG)
print("alpha_KMC avg:", avg_alpha_kmc)
print("alpha_phenom avg:", avg_alpha_phenom)
