# ============================================================================
# Figure 1  --  A: cells on the published cESFW embedding
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM is a continuous branching trajectory, not a discrete island.
#   Establishes the topology on which every later analysis is ordered.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A: cells on the published cESFW embedding
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 1. The human ICM is a continuous, bipotent developmental
#   transition.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/radley_original.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

OUT = "/Users/adelosan/cesfw_runs/radley_original/"
os.makedirs(OUT, exist_ok=True)

Sample_Info, Counts, cESFW_Genes, Embedding, Sample_Info_all, Embedding_all = ldy.load_embryo_data()
print(f"Embedding shape: {Embedding.shape} (this is the ORIGINAL published cESFW UMAP, not a reconstruction)")

petro = Sample_Info[Sample_Info.Dataset == "Petropoulos 2016"].copy()
STATES = ["Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi", "Hypo"]
mask = petro["Manual_Annotations"].isin(STATES)
idx_positions = [Sample_Info.index.get_loc(i) for i in petro.index[mask]]
sub_emb = Embedding[idx_positions]
sub_meta = petro.loc[mask].copy()
sub_meta["UMAP1"] = sub_emb[:, 0]
sub_meta["UMAP2"] = sub_emb[:, 1]

print(f"\nCells plotted: {len(sub_meta)}")
print(sub_meta["Manual_Annotations"].value_counts())

# Also grab TE/other neighboring states for context (so branch points have their full context)
context_states = ["Mural TE", "Polar TE", "Mid TE", "Early TE", "8-Cell", "Embryonic Disc", "cTB", "sTB"]
ctx_mask = petro["Manual_Annotations"].isin(context_states)
ctx_idx = [Sample_Info.index.get_loc(i) for i in petro.index[ctx_mask]]
ctx_emb = Embedding[ctx_idx]

fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(ctx_emb[:, 0], ctx_emb[:, 1], s=10, alpha=0.25, color="lightgrey", label="TE/other (context)")

colors = {"Morula": "#440154", "ICM/TE Branch": "#3b528b", "ICM": "#21918c",
          "Epi/Hypo Branch": "#5ec962", "preIm-Epi": "#fde725", "Hypo": "#e34a33"}
for s in STATES:
    d = sub_meta[sub_meta["Manual_Annotations"] == s]
    ax.scatter(d.UMAP1, d.UMAP2, s=35, alpha=0.85, color=colors[s], label=f"{s} (n={len(d)})",
               edgecolors="white", linewidths=0.3)

ax.set_xlabel("UMAP1 (original cESFW/Radley embedding)")
ax.set_ylabel("UMAP2 (original cESFW/Radley embedding)")
ax.set_title("Petropoulos 2016 cells on the ORIGINAL published cESFW branching topology\n(not a reconstruction -- this is Radley et al.'s own embedding)")
ax.legend(fontsize=9, loc="best")
plt.tight_layout()
plt.savefig(OUT + "radley_original_topology.png", dpi=150)
print(f"\nSaved: {OUT}radley_original_topology.png")
