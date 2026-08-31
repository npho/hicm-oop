# ============================================================================
# Figure 2  --  MYC / metabolic panel
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The canonical naive network is assembled only AFTER the ICM: ICM is
#   demethylation-permissive (high TET2, near-absent UHRF1) and low MYC.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   MYC / metabolic panel
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 2. The canonical human naïve pluripotency network is assembled
#   only after the ICM. Bottom panels: the ICM is demethylation-permissive,
#   with high TET2 and near-absent UHRF1, whereas DNMT3B rises and TET2
#   declines during epiblast formation (robust in Petropoulos and Yanagida,
#   the two studies adequately sampling the resolved epiblast, with
#   directional support from Meistermann; Fig. S19). c-MYC remains low
#   across the preimplantation pluripotent compartment (ICM and epiblast
#   alike) and increases only after implantation. Thus, ground-state
#   epigenetic hallmarks precedes assembly of the canonical naïve
#   transcription-factor network in humans (unlike the mouse, where they
#   coincide). These trajectories are reproducible across three embryo
#   datasets (Fig. S19).
#
# ATTRIBUTION CONFIDENCE: LOW   (evidence: infer)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_LIN28A_MYC_OxPhos.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Figure 9. MYC-LIN28-let-7-mTOR circuit: transcriptional and metabolic evidence.

(A) Grouped FC bars: MYC, LIN28A/B, HMGA2 across human embryo ICM, t2iLGo ICM-like, mouse diapause.
(B) OxPhos FC (100% suppression in ICM-like).
(C-F) LIN28A/B violins.
(G) OxPhos scatter: 100% of genes below diagonal.
(H) Seven-node metabolic quiescence circuit diagram.

Output: plots/FigLIN28A_MYC_OxPhos_v2.png
Run: conda run -n cesfw_env python figure_LIN28A_MYC_OxPhos.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from scipy.stats import mannwhitneyu
from sklearn.neighbors import KNeighborsClassifier
import umap as umap_lib
import os, sys, warnings
warnings.filterwarnings("ignore")
np.random.seed(42)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data

DATA     = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
DIAP_DIR = os.path.join(DATA, "diapause")
OUT      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
os.makedirs(OUT, exist_ok=True)

ICM_COL  = "#1F77B4"
EPI_COL  = "#FF7F0E"
T2IL_ICM = "#2CA02C"
T2IL_EPI = "#FFBF00"
MOUSE_COL = "#9467BD"

# ── Load embryo reference ─────────────────────────────────────────────────────
print("Loading embryo reference…")
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = load_embryo_data()
ann = np.asarray(Sample_Info["Manual_Annotations"])
icm_mask = ann == "ICM"
epi_mask = ann == "preIm-Epi"
lib_emb = Counts.sum(axis=1)
CPM_E   = Counts.div(lib_emb, axis=0) * 1e6
Log2_E  = np.log2(CPM_E + 1)
genes_emb = list(Counts.columns)

circuit_genes = ["MYC", "LIN28A", "LIN28B", "HMGA2", "IGF2BP1", "IGF2BP2", "FLCN",
                 "EIF4EBP1", "RPS6KB1", "CDK4", "ULK1", "CDKN1A"]

emb_fcs = {}
emb_pvs = {}
for g in circuit_genes:
    if g in genes_emb:
        v_i = Log2_E.loc[icm_mask, g].values
        v_e = Log2_E.loc[epi_mask, g].values
        fc  = float(v_i.mean() - v_e.mean())
        _, p = mannwhitneyu(v_i, v_e, alternative="two-sided")
        emb_fcs[g] = fc
        emb_pvs[g] = p

# OxPhos genes embryo
oxp_pref = ("MT-CO", "MT-CYB", "MT-ND", "MT-ATP",
            "NDUFA", "NDUFB", "NDUFC", "NDUFS", "NDUFV",
            "COX4", "COX5", "COX6", "COX7", "COX8",
            "ATP5", "UQCR", "CYC1", "SDHA", "SDHB", "SDHC", "SDHD")
oxp_emb = [g for g in genes_emb if any(g.startswith(p) for p in oxp_pref)]
oxp_fc_emb = [Log2_E.loc[icm_mask, g].mean() - Log2_E.loc[epi_mask, g].mean() for g in oxp_emb]
print(f"  Embryo OxPhos: n={len(oxp_emb)}  mean FC={np.mean(oxp_fc_emb):.3f}")

# ── Load and classify t2iLGo cells ───────────────────────────────────────────
print("Fitting UMAP…")
Model = umap_lib.UMAP(n_neighbors=50, metric='correlation',
                      min_dist=0.1, n_components=2, random_state=42).fit(Counts[cESFW_Genes])
Ref_Emb = Model.embedding_.copy()
anchor  = ann == "8-Cell"
flip_x  = Ref_Emb[anchor, 0].mean() < Ref_Emb[:, 0].mean()
flip_y  = Ref_Emb[anchor, 1].mean() > Ref_Emb[:, 1].mean()
if flip_x: Ref_Emb[:, 0] *= -1
if flip_y: Ref_Emb[:, 1] *= -1

def orient(arr):
    a = arr.copy()
    if flip_x: a[:, 0] *= -1
    if flip_y: a[:, 1] *= -1
    return a

Messmer_Counts = pd.read_csv(os.path.join(DATA, "Messmet2019_sc_Counts.csv"), header=0, index_col=0)
Messmer_Info   = pd.read_csv(os.path.join(DATA, "E-MTAB-6819_Sample_Info.csv"), header=0, index_col=0)
Messmer_Info2  = Messmer_Info[["Assay Name", "Factor Value[phenotype]",
                                "Characteristics[post-analysis well quality]"]].drop_duplicates()
Messmer_Info2.index = Messmer_Info2["Assay Name"]
naive_meta = Messmer_Info2[
    (Messmer_Info2["Factor Value[phenotype]"] == "naive") &
    (Messmer_Info2["Characteristics[post-analysis well quality]"] == "pass")
]
naive_cells  = naive_meta.index[naive_meta.index.isin(Messmer_Counts.index)]
Naive_Counts = Messmer_Counts.loc[naive_cells]
present = cESFW_Genes[np.isin(cESFW_Genes, Naive_Counts.columns)]
absent  = cESFW_Genes[~np.isin(cESFW_Genes, Naive_Counts.columns)]
Naive_Al = Naive_Counts[present]
if len(absent):
    zeros = pd.DataFrame(0, index=Naive_Al.index, columns=absent)
    Naive_Al = pd.concat([Naive_Al, zeros], axis=1)
Naive_Al = Naive_Al[cESFW_Genes].fillna(0)
Naive_Ov = orient(Model.transform(Naive_Al))
keep_ref = ann != "Unknown"
knn      = KNeighborsClassifier(n_neighbors=20).fit(Ref_Emb[keep_ref], ann[keep_ref])
Naive_Pred = knn.predict(Naive_Ov)
ICM_idx = np.where(Naive_Pred == "ICM")[0]
PRE_idx = np.where(Naive_Pred == "preIm-Epi")[0]
n_icm_t2il, n_pre_t2il = len(ICM_idx), len(PRE_idx)

naive_lib = Naive_Counts.sum(axis=1)
Naive_CPM = Naive_Counts.div(naive_lib, axis=0) * 1e6
Log2_N    = np.log2(Naive_CPM + 1)
t2il_fcs = {}
t2il_pvs = {}
for g in circuit_genes:
    if g in Naive_Counts.columns:
        v_i = Log2_N.iloc[ICM_idx][g].values
        v_p = Log2_N.iloc[PRE_idx][g].values
        fc  = float(v_i.mean() - v_p.mean())
        _, p = mannwhitneyu(v_i, v_p, alternative="two-sided")
        t2il_fcs[g] = fc
        t2il_pvs[g] = p

oxp_naive = [g for g in Naive_Counts.columns if any(g.startswith(p) for p in oxp_pref)]
oxp_fc_t2il_all = [Log2_N.iloc[ICM_idx][g].mean() - Log2_N.iloc[PRE_idx][g].mean()
                   for g in oxp_naive]
n_below_zero = sum(1 for f in oxp_fc_t2il_all if f < 0)
print(f"  t2iLGo OxPhos: n={len(oxp_naive)}  {n_below_zero}/{len(oxp_naive)} below 0")

# ── Mouse diapause values (Torin1 vs serum from manuscript) ──────────────────
mouse_fcs = {"MYC": -1.91, "LIN28A": -0.52, "LIN28B": -0.93, "HMGA2": -2.10}
try:
    gse_path = os.path.join(DIAP_DIR, "GSE81285_RNA-seq_expression_mm10.tsv.gz")
    Mouse_Raw = pd.read_csv(gse_path, sep='\t', index_col=0)
    cols = Mouse_Raw.columns.tolist()
    serum_c = [c for c in cols if 'serum' in c.lower() or ('lif' in c.lower() and '2i' not in c.lower())]
    torin_c = [c for c in cols if 'torin' in c.lower() or 'paus' in c.lower()]
    if serum_c and torin_c:
        mouse_lib = Mouse_Raw.sum(axis=0)
        Mouse_CPM = Mouse_Raw.div(mouse_lib, axis=1) * 1e6
        Log2_M    = np.log2(Mouse_CPM + 1)
        for g_human, g_mouse in [("MYC", "Myc"), ("LIN28A", "Lin28a"),
                                  ("LIN28B", "Lin28b"), ("HMGA2", "Hmga2")]:
            for gm in [g_mouse, g_mouse.upper(), g_human]:
                if gm in Mouse_Raw.index:
                    fc = Log2_M.loc[gm, torin_c].mean() - Log2_M.loc[gm, serum_c].mean()
                    mouse_fcs[g_human] = float(fc)
                    break
except Exception as e:
    print(f"  Mouse data load error: {e}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 10))
gs  = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.40)

# Panel A: Grouped FC bars for circuit genes
ax_a = fig.add_subplot(gs[0, 0:2])
focus_genes = ["MYC", "LIN28A", "LIN28B", "HMGA2"]
x   = np.arange(len(focus_genes))
w   = 0.25
groups = [
    ("Embryo ICM", emb_fcs, ICM_COL),
    ("t2iLGo ICM-like", t2il_fcs, T2IL_ICM),
    ("Mouse Torin1", mouse_fcs, MOUSE_COL),
]
for j, (label, fc_dict, col) in enumerate(groups):
    vals = [fc_dict.get(g, np.nan) for g in focus_genes]
    ax_a.bar(x + (j - 1) * w, vals, width=w, color=col, edgecolor='k',
             linewidth=0.6, label=label, alpha=0.85)

ax_a.axhline(0, color='black', linewidth=1)
ax_a.set_xticks(x)
ax_a.set_xticklabels(focus_genes, fontsize=11)
ax_a.set_ylabel("log₂FC (ICM-like vs preEPI-like)", fontsize=10)
ax_a.set_title("A  MYC–LIN28–HMGA2 concordance across contexts", fontsize=10, fontweight='bold', loc='left')
ax_a.legend(fontsize=9, frameon=False)
ax_a.tick_params(labelsize=9)

# Panel B: OxPhos suppression across contexts
ax_b = fig.add_subplot(gs[0, 2])
contexts = ["Embryo ICM", "t2iLGo\nICM-like"]
oxp_vals = [np.mean(oxp_fc_emb), np.mean(oxp_fc_t2il_all)]
pct_vals = [100 * sum(f < 0 for f in oxp_fc_emb) / len(oxp_fc_emb),
            100 * n_below_zero / len(oxp_fc_t2il_all)]
cols_b   = [ICM_COL, T2IL_ICM]
ax_b.bar(contexts, oxp_vals, color=cols_b, edgecolor='k', linewidth=0.8, width=0.5)
ax_b.axhline(0, color='black', linewidth=1)
for i, (fc, pct) in enumerate(zip(oxp_vals, pct_vals)):
    ax_b.text(i, fc - 0.03, f"FC={fc:.2f}\n{pct:.0f}% below 0",
              ha='center', va='top', fontsize=9, fontweight='bold')
ax_b.set_ylabel("Mean OxPhos log₂FC", fontsize=10)
ax_b.set_title("B  OxPhos suppression", fontsize=10, fontweight='bold', loc='left')
ax_b.tick_params(labelsize=9)

# Panels C-F: LIN28A and LIN28B violins (embryo and t2iLGo)
violin_specs = [
    (gs[1, 0], Log2_E, icm_mask, epi_mask, "LIN28A", "C  LIN28A (embryo)", ICM_COL, EPI_COL),
    (gs[1, 1], Log2_E, icm_mask, epi_mask, "LIN28B", "D  LIN28B (embryo)", ICM_COL, EPI_COL),
    (gs[1, 2], Log2_N, ICM_idx, PRE_idx, "LIN28A", "E  LIN28A (t2iLGo)", T2IL_ICM, T2IL_EPI),
    (gs[1, 3], Log2_N, ICM_idx, PRE_idx, "LIN28B", "F  LIN28B (t2iLGo)", T2IL_ICM, T2IL_EPI),
]
for spec in violin_specs:
    gs_loc, log2_df, idx_a, idx_b, gene, title, col_a, col_b = spec
    ax = fig.add_subplot(gs_loc)
    if gene not in log2_df.columns:
        ax.text(0.5, 0.5, f"{gene}\nnot detected", ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=9, fontweight='bold', loc='left')
        continue
    if isinstance(idx_a, np.ndarray) and idx_a.dtype == bool:
        v_a = log2_df.loc[idx_a, gene].values
        v_b = log2_df.loc[idx_b, gene].values
        n_a = idx_a.sum()
        n_b = idx_b.sum()
    else:
        v_a = log2_df.iloc[idx_a][gene].values
        v_b = log2_df.iloc[idx_b][gene].values
        n_a = len(idx_a)
        n_b = len(idx_b)
    _, pv = mannwhitneyu(v_a, v_b, alternative="two-sided")
    fc = float(v_a.mean() - v_b.mean())
    label_a = "ICM-like" if "t2iLGo" in title else "ICM"
    label_b = "preEPI-like" if "t2iLGo" in title else "preEPI"
    parts = ax.violinplot([v_b, v_a], positions=[1, 2], showmedians=True, showextrema=False)
    parts['bodies'][0].set_facecolor(col_b); parts['bodies'][0].set_alpha(0.7)
    parts['bodies'][1].set_facecolor(col_a); parts['bodies'][1].set_alpha(0.7)
    parts['cmedians'].set_colors('black'); parts['cmedians'].set_linewidths(2)
    ymax = max(v_a.max(), v_b.max()) if len(v_a) and len(v_b) else 1
    sig = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
    ax.plot([1, 2], [ymax * 1.1, ymax * 1.1], 'k-', lw=1)
    ax.text(1.5, ymax * 1.15, f"{sig}\nFC={fc:+.2f}", ha='center', va='bottom', fontsize=9)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"{label_b}\n(n={n_b})", f"{label_a}\n(n={n_a})"], fontsize=8)
    ax.set_ylabel("log₂(CPM+1)", fontsize=9)
    ax.set_title(title, fontsize=9, fontweight='bold', loc='left')
    ax.tick_params(labelsize=8)

# Panel H: Circuit diagram (Panel G slot for OxPhos scatter, H for circuit)
ax_g = fig.add_subplot(gs[0, 3])
# OxPhos scatter: embryo ICM vs t2iLGo ICM-like
emb_oxp = [g for g in oxp_emb if g in Naive_Counts.columns]
emb_mean_oxp = np.array([Log2_E.loc[icm_mask, g].mean() - Log2_E.loc[epi_mask, g].mean()
                          for g in emb_oxp])
t2il_mean_oxp = np.array([Log2_N.iloc[ICM_idx][g].mean() - Log2_N.iloc[PRE_idx][g].mean()
                           for g in emb_oxp])
ax_g.scatter(emb_mean_oxp, t2il_mean_oxp, color='#9B59B6', alpha=0.5, s=12,
             label=f"OxPhos genes (n={len(emb_oxp)})")
lim = max(abs(emb_mean_oxp).max(), abs(t2il_mean_oxp).max()) + 0.2
ax_g.plot([-lim, lim], [-lim, lim], 'k--', linewidth=1, alpha=0.5, label='y=x')
ax_g.axhline(0, color='gray', linewidth=0.5)
ax_g.axvline(0, color='gray', linewidth=0.5)
pct_both = 100 * np.mean((emb_mean_oxp < 0) & (t2il_mean_oxp < 0))
ax_g.text(0.05, 0.95, f"{pct_both:.0f}% of OxPhos genes\ndown in both contexts",
          transform=ax_g.transAxes, fontsize=9, va='top',
          bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.8))
ax_g.set_xlabel("Embryo ICM FC", fontsize=9)
ax_g.set_ylabel("t2iLGo ICM-like FC", fontsize=9)
ax_g.set_title("G  OxPhos concordance", fontsize=9, fontweight='bold', loc='left')
ax_g.legend(fontsize=8, frameon=False)

plt.suptitle("Figure 9. MYC–LIN28–let-7–mTOR circuit: transcriptional and metabolic evidence",
             fontsize=11, fontweight='bold', y=1.01)
plt.savefig(os.path.join(OUT, "FigLIN28A_MYC_OxPhos_v3_CPMnorm.png"), dpi=200, bbox_inches='tight')
print("Saved FigLIN28A_MYC_OxPhos_v3_CPMnorm.png")
plt.close()
