# ============================================================================
# Figure 1  --  PAGA branch-hub / betweenness panel
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM is a continuous branching trajectory, not a discrete island.
#   Establishes the topology on which every later analysis is ordered.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   PAGA branch-hub / betweenness panel
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 1. The human ICM is a continuous, bipotent developmental
#   transition.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/Fig1_PAGA_branch_hub.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Fig 1 — PAGA connectivity on the cESFW human-embryo atlas.

Backs the manuscript caption ("Partition-based graph abstraction on the stage
annotations; ... the ICM as the branch hub"). PAGA partitions = the cESFW-atlas
Manual_Annotations; connectivity computed from a correlation-metric KNN graph in
the cESFW feature-gene space (project convention). Node positions initialised
from the published cESFW UMAP embedding so the graph is spatially faithful.

Run: conda run -n cesfw_env python Fig1_PAGA_branch_hub.py
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")
os.makedirs(OUT, exist_ok=True)

# Developmental ordering of the cESFW atlas annotations (for stable colouring/layout)
ORDER = ["8-Cell", "Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch",
         "preIm-Epi", "Embryonic Disc", "Hypo", "ExE-Mech",
         "Early TE", "Mid TE", "Polar TE", "Mural TE", "cTB", "sTB"]

Sample_Info, Counts, cESFW_Genes, Embedding, SI_full, Emb_full = load_embryo_data()
ann = Sample_Info["Manual_Annotations"].astype(str)
keep = (ann != "Unknown").values
feat = [g for g in cESFW_Genes if g in Counts.columns]

lib = Counts.sum(1); lib[lib == 0] = 1
log2cpm = np.log2(Counts.div(lib, axis=0) * 1e6 + 1)

adata = sc.AnnData(X=log2cpm.loc[keep, feat].values)
cats = [c for c in ORDER if c in set(ann[keep])]
adata.obs["annotation"] = pd.Categorical(ann[keep].values, categories=cats, ordered=True)
adata.obsm["X_umap"] = Embedding[keep]

print(f"PAGA on {adata.n_obs} annotated cells x {adata.n_vars} cESFW genes; "
      f"{len(cats)} groups")

sc.pp.neighbors(adata, n_neighbors=15, use_rep="X", metric="correlation")
sc.tl.paga(adata, groups="annotation")

# ── connectivity table ──────────────────────────────────────────────────────
conn = pd.DataFrame(adata.uns["paga"]["connectivities"].toarray(),
                    index=cats, columns=cats)
conn.to_csv(os.path.join(OUT, "PAGA_connectivities.csv"))

# degree-weighted "hub" score for each node (sum of connectivities)
hub = conn.sum(1).sort_values(ascending=False)
print("\nPAGA connectivity degree (sum of edge weights):")
print(hub.round(3).to_string())
hub.to_csv(os.path.join(OUT, "PAGA_hub_degree.csv"))

# ── betweenness centrality: is the ICM the BRIDGE/internode? ─────────────────
# (degree measures how connected; betweenness measures how much a node lies on
#  the shortest paths between all others — the correct test for "internode".)
import networkx as nx
G = nx.Graph()
G.add_nodes_from(cats)
thr = 0.05
for i, a in enumerate(cats):
    for j in range(i + 1, len(cats)):
        w = conn.iloc[i, j]
        if w > thr:
            # distance = inverse connectivity (strong edge -> short path)
            G.add_edge(a, cats[j], weight=w, distance=1.0 / w)
btw = pd.Series(nx.betweenness_centrality(G, weight="distance", normalized=True))
btw = btw.sort_values(ascending=False)
print("\nBetweenness centrality (bridge score; higher = more of an internode):")
print(btw.round(3).to_string())
btw.to_csv(os.path.join(OUT, "PAGA_betweenness.csv"))
print(f"\n>> [ALL 15 states] Top internode node: {btw.index[0]} "
      f"(betweenness={btw.iloc[0]:.3f}); ICM rank = "
      f"{list(btw.index).index('ICM') + 1}/{len(btw)}")

# ── lineage-restricted betweenness: the claim is about the EMBRYONIC (ICM-derived)
#    lineage, so recompute on that subgraph only (drop TE terminals + ExE-Mech). ──
EMB_LINEAGE = ["8-Cell", "Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch",
               "preIm-Epi", "Embryonic Disc", "Hypo"]
emb_nodes = [c for c in EMB_LINEAGE if c in cats]
Gsub = G.subgraph(emb_nodes).copy()
btw_sub = pd.Series(nx.betweenness_centrality(Gsub, weight="distance",
                                              normalized=True)).sort_values(ascending=False)
print("\nBetweenness centrality — EMBRYONIC LINEAGE only "
      "(8-Cell..Hypo, TE/ExE excluded):")
print(btw_sub.round(3).to_string())
btw_sub.to_csv(os.path.join(OUT, "PAGA_betweenness_embryonic_lineage.csv"))
print(f"\n>> [Embryonic lineage] Top internode node: {btw_sub.index[0]} "
      f"(betweenness={btw_sub.iloc[0]:.3f}); ICM rank = "
      f"{list(btw_sub.index).index('ICM') + 1}/{len(btw_sub)}")

# ── node positions from UMAP centroids (spatially faithful layout) ───────────
cent = {}
for c in cats:
    m = (adata.obs["annotation"] == c).values
    cent[c] = adata.obsm["X_umap"][m].mean(0)
pos = np.array([cent[c] for c in cats])
adata.uns["paga"]["pos"] = pos

# ── figure ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))

sc.pl.paga(adata, threshold=0.05, node_size_scale=2.5, edge_width_scale=1.2,
           fontsize=8, ax=axes[0], show=False, frameon=True, pos=pos,
           title="PAGA connectivity (cESFW atlas, UMAP layout)\ngroups = stage annotations")

# UMAP coloured by annotation for context
sc.pl.umap(adata, color="annotation", ax=axes[1], show=False, size=18,
           legend_loc="right margin", title="cESFW UMAP — atlas annotations",
           frameon=True)

fig.tight_layout()
out = os.path.join(OUT, "Fig1_PAGA_branch_hub.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print("\nSaved:", out)
