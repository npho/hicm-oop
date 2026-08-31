# ============================================================================
# Figure S4  --  mouse E3.5->E4.5 sterol/mTOR delta-scores
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mTOR-sterol programme also rises across normal mouse E3.5 -> E4.5
#   (cross-species support that the axis is not a human artefact).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   mouse E3.5->E4.5 sterol/mTOR delta-scores
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S4. The mTOR–sterol programme rises across normal mouse
#   E3.5 → E4.5 (secondary, cross-species support). Sterol/SREBF (top) and
#   mTORC1 (bottom) delta-scores from the E3.5 ICM (common origin) into the
#   E4.5 epiblast and primitive-endoderm (PrE) fates, in Mohammed 2017
#   (Smart-seq2) and Nowotschin 2019 (10x). The conserved sterol increase is
#   carried by the PrE (hypoblast) arm; the epiblast-specific sterol rise is
#   weak and not robust across both datasets, so this supports the overall
#   ICM maturation rather than the epiblast lineage specifically.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_mouse_E35_E45.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Figure: sterol/SREBF + mTORC1 programs rise across normal mouse E3.5->E4.5,
two independent single-cell datasets, per lineage."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "plots", "revision_results", "mouse_E35_E45")
moh = pd.read_csv(os.path.join(OUT, "Mohammed_cell_scores.csv"), index_col=0)
now = pd.read_csv(os.path.join(OUT, "Nowotschin_cell_scores.csv"), index_col=0)

DATASETS = [("Mohammed 2017 (Smart-seq2)", moh), ("Nowotschin 2019 (10x)", now)]
PROGS = [("sterol", "Sterol / SREBF program"), ("mtorc1", "mTORC1 program")]
# E3.5 ICM = common origin (denominator); E4.5 resolved into the two fates
GROUPS = [("E3.5\nICM", "E3.5", None, "#8FB8DE"),
          ("E4.5\nEpiblast", "E4.5", "Epiblast", "#D64550"),
          ("E4.5\nPrE", "E4.5", "PrE", "#E8A33A")]

fig, axes = plt.subplots(2, 2, figsize=(11.5, 9))
for r, (prog, ptitle) in enumerate(PROGS):
    for c, (dname, df) in enumerate(DATASETS):
        ax = axes[r, c]
        d = df[df.lineage != "TE"]                       # ICM-derived only
        positions, ticklab = [], []
        for x, (glab, stg, lin, col) in enumerate(GROUPS):
            m = (d.stage == stg) if stg == "E3.5" else ((d.stage == stg) & (d.lineage == lin))
            v = d[m][prog].values
            bp = ax.boxplot(v, positions=[x], widths=0.6, patch_artist=True,
                            showfliers=False, medianprops=dict(color="black", lw=1.6))
            bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.85)
            jit = np.random.default_rng(0).normal(0, 0.06, len(v))
            ax.scatter(np.full(len(v), x) + jit, v, s=5, c="black", alpha=0.28, linewidths=0)
            positions.append(x); ticklab.append(f"{glab}\nn={len(v)}")
        ax.set_xticks(positions); ax.set_xticklabels(ticklab, fontsize=8.5)
        ax.axhline(0, color="grey", lw=0.6, ls=":")
        ax.axvline(0.5, color="grey", lw=0.5, ls="--", alpha=0.6)  # origin | fates
        if c == 0: ax.set_ylabel(f"{ptitle}\n(delta-score)", fontsize=11)
        if r == 0: ax.set_title(dname, fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)

from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor="#8FB8DE", label="E3.5 ICM (common origin)"),
                    Patch(facecolor="#D64550", label="E4.5 epiblast"),
                    Patch(facecolor="#E8A33A", label="E4.5 PrE")],
           loc="upper center", ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 1.02))
fig.suptitle("Sterol/SREBF and mTORC1 programs rise across normal mouse E3.5 → E4.5\n"
             "(two independent single-cell datasets, per lineage)", fontsize=13, y=1.07)
fig.tight_layout()
out = os.path.join(OUT, "Mouse_E35_E45_sterol_mTORC1.png")
fig.savefig(out, dpi=250, bbox_inches="tight")
print("Saved:", out)
