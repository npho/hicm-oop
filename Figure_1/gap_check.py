# ============================================================================
# Figure 1  --  continuity check across the ICM/TE-Epi/Hypo axis (no gap)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM is a continuous branching trajectory, not a discrete island.
#   Establishes the topology on which every later analysis is ordered.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   continuity check across the ICM/TE-Epi/Hypo axis (no gap)
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 1. The human ICM is a continuous, bipotent developmental
#   transition.
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: infer)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/gap_check.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

OUT = "/Users/adelosan/cesfw_runs/proper_trajectory/"
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()

# Use the ORIGINAL published cESFW embedding (Radley), not a recomputation
petro = Sample_Info[Sample_Info.Dataset == "Petropoulos 2016"].copy()
STATES = ["Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
mask = petro["Manual_Annotations"].isin(STATES)
idx_pos = [Sample_Info.index.get_loc(i) for i in petro.index[mask]]
sub = petro.loc[mask].copy()
sub["U1"] = Embedding[idx_pos, 0]
sub["U2"] = Embedding[idx_pos, 1]

colors = {"Morula":"#3b7dd8","ICM/TE Branch":"#f59e0b","ICM":"#16a34a",
          "Epi/Hypo Branch":"#dc2626","preIm-Epi":"#9333ea","Hypo":"#8b5a2b"}

# ── Define the axis connecting the two branch-point centroids ──────────────
# ICM/TE Branch centroid -> Epi/Hypo Branch centroid, project all cells onto it
c_icmte = sub.loc[sub.Manual_Annotations=="ICM/TE Branch", ["U1","U2"]].mean().values
c_epihypo = sub.loc[sub.Manual_Annotations=="Epi/Hypo Branch", ["U1","U2"]].mean().values
axis = c_epihypo - c_icmte
axis = axis / np.linalg.norm(axis)
sub["proj"] = (sub[["U1","U2"]].values - c_icmte) @ axis  # 0 at ICM/TE centroid, |axis-length| at Epi/Hypo

d_branches = np.linalg.norm(c_epihypo - c_icmte)
print(f"Distance between ICM/TE Branch and Epi/Hypo Branch centroids (embedding units): {d_branches:.2f}")
print("\nProjection position of each stage centroid along the ICM/TE -> Epi/Hypo axis:")
for s in STATES:
    v = sub.loc[sub.Manual_Annotations==s, "proj"]
    print(f"  {s:16s}: mean={v.mean():+.2f}  (n={len(v)})")

# ── Count cells in the interval strictly BETWEEN the two branch centroids ──
between = sub[(sub["proj"] > 0.15*d_branches) & (sub["proj"] < 0.85*d_branches)]
print(f"\nCells in the middle 70% of the interval between the two branch points: {len(between)}")
print("  Their stage labels:", dict(between.Manual_Annotations.value_counts()))

# ── Histogram of density along that axis: is there a trough (gap) or is it filled? ──
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Left: the real embedding, cells colored by stage, with the two branch centroids marked
ax = axes[0]
# context (grey) - all other petro cells
ctx = petro[~mask]
ctx_pos = [Sample_Info.index.get_loc(i) for i in ctx.index]
ax.scatter(Embedding[ctx_pos,0], Embedding[ctx_pos,1], s=8, c="#e5e5e5", alpha=0.4)
for s in STATES:
    d = sub[sub.Manual_Annotations==s]
    ax.scatter(d.U1, d.U2, s=28, color=colors[s], alpha=0.85, edgecolors="white", linewidths=0.3, label=s)
ax.scatter(*c_icmte, s=300, marker="X", color="black", zorder=10)
ax.scatter(*c_epihypo, s=300, marker="X", color="black", zorder=10)
ax.plot([c_icmte[0], c_epihypo[0]], [c_icmte[1], c_epihypo[1]], "k--", lw=1.5, zorder=9)
ax.legend(fontsize=8, loc="best")
ax.set_title("Original published cESFW embedding\n(X = branch-point centroids)")
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")

# Right: density histogram along the connecting axis
ax = axes[1]
bins = np.linspace(sub["proj"].min(), sub["proj"].max(), 40)
for s in STATES:
    d = sub[sub.Manual_Annotations==s]
    ax.hist(d["proj"], bins=bins, color=colors[s], alpha=0.7, label=s, stacked=False)
ax.axvline(0, color="orange", ls="--", lw=2, label="ICM/TE centroid")
ax.axvline(d_branches, color="red", ls="--", lw=2, label="Epi/Hypo centroid")
ax.set_xlabel("Position along ICM/TE Branch -> Epi/Hypo Branch axis")
ax.set_ylabel("Cell count")
ax.set_title("Is the region between the two branches filled, or a gap?")
ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig(OUT + "gap_check.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: {OUT}gap_check.png")

# ── Formal: nearest-neighbor gap test. Sort cells by projection, look at spacing ──
proj_sorted = np.sort(sub["proj"].values)
gaps = np.diff(proj_sorted)
print(f"\nLargest gaps in cell spacing along the axis (should be small/uniform if continuous):")
big_gap_idx = np.argsort(gaps)[-5:]
for gi in sorted(big_gap_idx):
    print(f"  gap of {gaps[gi]:.2f} between proj={proj_sorted[gi]:.2f} and proj={proj_sorted[gi+1]:.2f}")
print(f"  median gap={np.median(gaps):.3f}, mean gap={np.mean(gaps):.3f}, max gap={gaps.max():.3f}")
print(f"  max/median gap ratio: {gaps.max()/np.median(gaps):.1f} (high ratio = a real discontinuity)")
