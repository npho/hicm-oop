# ============================================================================
# Figure S1  --  stress/apoptosis/UPR scores + QC metrics by lineage
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Branch/ICM cells are not stressed or apoptotic 'reject' cells - rules
#   out the main alternative explanation for a transcriptionally
#   intermediate population.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   stress/apoptosis/UPR scores + QC metrics by lineage
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S1. cESFW branch/ICM cells show no distinct stress/apoptosis
#   “reject” signature. Apoptosis and stress (unfolded-protein, heat-shock
#   and oxidative) programme delta-scores, percentage mitochondrial reads,
#   and genes detected per cell across human embryo lineages in
#   developmental order; the questioned cESFW groups (ICM/TE Branch, ICM,
#   Epi/Hypo Branch; red) versus other lineages (grey). Apoptosis scores in
#   the questioned groups lie within the range of bona fide stages (8-cell,
#   morula and embryonic disc are comparable or higher); stress programmes
#   are not elevated (if anything lower); the modestly lower complexity of
#   the ICM and Epi/Hypo Branch is consistent with the low transcriptional
#   output of the ICM and largely reflects library size. The cells therefore
#   lack the coordinated high-stress, high-apoptosis, low-quality profile
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_reject_cells.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Figure: are cESFW branch/ICM cells 'REJECT' (stressed/apoptotic/low-quality)?
Stress/apoptosis programs + QC across lineages (developmental order); questioned
groups highlighted."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")
df = pd.read_csv(os.path.join(OUT, "reject_cell_scores.csv"), index_col=0)

ORDER = ["8-Cell","Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo",
         "Embryonic Disc","Early TE","Mid TE","Polar TE","Mural TE","cTB","sTB","ExE-Mech"]
ORDER = [o for o in ORDER if o in df.ann.unique()]
QUESTIONED = {"ICM/TE Branch","ICM","Epi/Hypo Branch"}
PANELS = [("apoptosis","Apoptosis program (delta-score)"),
          ("stress","Stress / UPR / oxidative (delta-score)"),
          ("pct_mito","% mitochondrial reads"),
          ("n_genes","Genes detected per cell")]

fig, axes = plt.subplots(2, 2, figsize=(15, 10))
for ax, (metric, title) in zip(axes.ravel(), PANELS):
    for x, lin in enumerate(ORDER):
        v = df[df.ann == lin][metric].values
        col = "#D64550" if lin in QUESTIONED else "#9DB4C0"
        bp = ax.boxplot(v, positions=[x], widths=0.64, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.5))
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.9)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=55, ha="right", fontsize=8.5)
    ax.set_title(title, fontsize=12)
    ax.spines[["top","right"]].set_visible(False)

from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor="#D64550", label="Questioned cESFW branch/ICM groups"),
                    Patch(facecolor="#9DB4C0", label="Other lineages")],
           loc="upper center", ncol=2, frameon=False, fontsize=11, bbox_to_anchor=(0.5, 1.015))
fig.suptitle("cESFW branch/ICM cells show no distinct stress/apoptosis 'reject' signature",
             fontsize=13.5, y=1.05)
fig.tight_layout()
out = os.path.join(OUT, "Reject_cells_stress_QC.png")
fig.savefig(out, dpi=250, bbox_inches="tight")
print("Saved:", out)
