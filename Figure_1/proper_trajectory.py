# ============================================================================
# Figure 1  --  B,C: PAGA graph + diffusion pseudotime; writes the *_with_dpt.csv used downstream
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM is a continuous branching trajectory, not a discrete island.
#   Establishes the topology on which every later analysis is ordered.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   B,C: PAGA graph + diffusion pseudotime; writes the *_with_dpt.csv used
#   downstream
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
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/proper_trajectory.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Proper, branch-aware trajectory analysis using scanpy's PAGA + diffusion
pseudotime (DPT) on the SAME gene set Radley et al.'s cESFW topology was
built on (Saved_cESFW_Genes.npy) -- not a from-scratch PCA on a hand-picked
subset. This natively handles the two-branch topology (Morula -> ICM/TE
Branch -> ICM -> Epi/Hypo Branch -> {preEPI, Hypo}) instead of flattening
it onto one linear axis.

Then: re-test mTORC1 signaling and EPI/Hypo lineage bias as a function of
DPT pseudotime (embryo-level, circularity-free pathway scores as before),
properly this time.

Run: conda run -n cesfw_env python proper_trajectory.py
"""
import numpy as np, pandas as pd, sys, os
import scanpy as sc
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from scipy.stats import spearmanr
import statsmodels.api as sm

OUT = "/Users/adelosan/cesfw_runs/proper_trajectory/"
os.makedirs(OUT, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
petro = Sample_Info[Sample_Info.Dataset == "Petropoulos 2016"].copy()
STATES = ["Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi", "Hypo"]
sub = petro[petro["Manual_Annotations"].isin(STATES)].copy()
counts_sub = Counts.loc[sub.index]

feat_genes = [g for g in cESFW_Genes if g in counts_sub.columns]
print(f"Using the ORIGINAL {len(feat_genes)} cESFW/Radley genes for trajectory structure "
      f"(same gene space as the published topology, not a hand-picked subset)")

# ── Build AnnData, standard scanpy preprocessing on the cESFW gene space ──
adata = ad.AnnData(X=counts_sub[feat_genes].values.astype(np.float32))
adata.obs = sub.reset_index().rename(columns={"index": "cell_id"})
adata.var_names = feat_genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, random_state=0)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(adata, random_state=0)

# ── PAGA: confirm branch connectivity between the labeled stages ──────────
adata.obs["Manual_Annotations"] = pd.Categorical(adata.obs["Manual_Annotations"], categories=STATES)
sc.tl.paga(adata, groups="Manual_Annotations")
print("\nPAGA connectivity matrix (rows/cols = stages, values = connection strength 0-1):")
paga_df = pd.DataFrame(adata.uns["paga"]["connectivities"].toarray(),
                        index=STATES, columns=STATES)
print(paga_df.round(2).to_string())
paga_df.to_csv(OUT + "paga_connectivity.csv")

fig, ax = plt.subplots(figsize=(7, 6))
sc.pl.paga(adata, ax=ax, show=False, node_size_scale=2, fontsize=8)
plt.savefig(OUT + "paga_graph.png", dpi=150, bbox_inches="tight")
print(f"Saved: {OUT}paga_graph.png")

# ── Diffusion pseudotime, rooted at a Morula cell ──────────────────────────
sc.tl.diffmap(adata)
morula_idx = np.where(adata.obs["Manual_Annotations"] == "Morula")[0]
# root = the morula cell most extreme along diffmap component 1
root = morula_idx[np.argmax(np.abs(adata.obsm["X_diffmap"][morula_idx, 1]))]
adata.uns["iroot"] = int(root)
sc.tl.dpt(adata)
sub["dpt_pseudotime"] = adata.obs["dpt_pseudotime"].values

print("\nDPT pseudotime by stage (should increase Morula -> ... -> Hypo/preEPI):")
print(sub.groupby("Manual_Annotations")["dpt_pseudotime"].agg(["mean", "std", "count"]).reindex(STATES))

# ── UMAP plot colored by stage and by pseudotime ───────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="Manual_Annotations", ax=axes[0], show=False, title="Scanpy UMAP (cESFW genes) -- stage")
sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[1], show=False, title="DPT pseudotime (rooted at Morula)")
plt.tight_layout()
plt.savefig(OUT + "umap_stage_pseudotime.png", dpi=150)
print(f"Saved: {OUT}umap_stage_pseudotime.png")

sub.to_csv(OUT + "petropoulos_with_dpt.csv")
adata.write(OUT + "adata.h5ad")
print(f"\nSaved AnnData + cell-level pseudotime table to {OUT}")
print("\nDone with structure inference. Run part 2 for pathway re-analysis.")
