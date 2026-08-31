# ============================================================================
# Figure S6  --  cynomolgus analysis
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ESRRB ICM-stage rewiring is conserved in a second primate
#   (cynomolgus).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   cynomolgus analysis
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S6. ρ=-0.81ρ=+0.77The ESRRB ICM-stage rewiring is conserved in a
#   second primate (cynomolgus). cESFW applied de novo to cynomolgus pre-
#   implantation cells (Nakamura et al., 2016; GSE74767), the same framework
#   used for the human atlas. (A) cESFW recovers a branch topology: an
#   early/ICM–TE centre resolving to epiblast and hypoblast. (B) Diffusion
#   pseudotime (rooted in the early/ICM compartment) flows toward both
#   fates. (C) ESRRB is high in the early/ICM compartment and low in the
#   committed epiblast; along the epiblast lineage ESRRB declines (Spearman,
#   n=45) while the naïve network rises (), as in human. Caveats: this is a
#   de novo marker-annotated re-analysis; ESRRB is high in the early
#   trophectoderm too (so it marks the early compartment, not the ICM
#   exclusively, as in human); the ICM as a distinct populated cluster is
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_nhp_cynomolgus.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
NHP Preimplantation Embryo UMAP Analysis — Cynomolgus Macaque
Dataset: GSE74767, Nakamura et al. 2017 Scientific Data
420 cells total; using E06-E09 (193 preimplantation cells)
Gene symbols: macFas5 (largely matching human HGNC)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import gzip
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from umap import UMAP
import warnings
warnings.filterwarnings('ignore')

# ── 1. Load data ──────────────────────────────────────────────────────────────
print("Loading cynomolgus data...")
DATA = "/Users/adelosan/cESFW_analysis/data/nhp/cynomolgus/GSE74767_SC3seq_Cy_ProcessedData.txt.gz"

with gzip.open(DATA, 'rt') as f:
    raw = pd.read_csv(f, sep='\t', index_col=None)

# columns: macFas5_entrez_id | macFas5_gene_symbol | cell1 | cell2 ...
gene_symbols = raw.iloc[:, 1].values
cell_cols = raw.columns[2:].tolist()
expr_all = raw.iloc[:, 2:].values.astype(float)
expr_df = pd.DataFrame(expr_all, index=gene_symbols, columns=cell_cols)

# Remove duplicate gene symbols (keep max)
expr_df = expr_df.groupby(expr_df.index).max()

print(f"Full dataset: {expr_df.shape[0]} genes × {expr_df.shape[1]} cells")

# ── 2. Select E06-E09 preimplantation cells ───────────────────────────────────
preimpl_stages = ['E06', 'E07', 'E08', 'E09']
preimpl_cells = [c for c in cell_cols if c[:3] in preimpl_stages]
expr = expr_df[preimpl_cells].copy()
stage_map = {c: c[:3] for c in preimpl_cells}
stages = pd.Series([stage_map[c] for c in preimpl_cells], index=preimpl_cells, name='stage')

print(f"Preimplantation cells: {expr.shape[1]}")
print(stages.value_counts().sort_index())

# ── 3. Log2(RPKM+1) transform ─────────────────────────────────────────────────
log_expr = np.log2(expr + 1)

# ── 4. Filter and select highly variable genes ────────────────────────────────
# Remove genes expressed in <3 cells
min_cells = 3
expressed = (expr > 0).sum(axis=1)
log_expr = log_expr[expressed >= min_cells]
print(f"After filtering low-expressed genes: {log_expr.shape[0]} genes")

# Select top 2000 HVGs by dispersion (variance / mean)
gene_mean = log_expr.mean(axis=1)
gene_var = log_expr.var(axis=1)
gene_cv2 = gene_var / (gene_mean + 1e-6)
hvgs = gene_cv2.nlargest(2000).index.tolist()
X = log_expr.loc[hvgs].T.values  # cells × genes

print(f"Using {len(hvgs)} HVGs for UMAP")

# ── 5. Build UMAP ─────────────────────────────────────────────────────────────
print("Running UMAP...")
reducer = UMAP(n_components=2, n_neighbors=15, min_dist=0.3,
               metric='euclidean', random_state=42)
embedding = reducer.fit_transform(X)
umap_df = pd.DataFrame(embedding, index=preimpl_cells, columns=['UMAP1', 'UMAP2'])
umap_df['stage'] = stages.values

# ── 6. Marker-based lineage annotation ───────────────────────────────────────
def score_gene(gene, cells):
    if gene in log_expr.index:
        return log_expr.loc[gene, cells].values
    return np.zeros(len(cells))

nanog   = score_gene('NANOG', preimpl_cells)
gata6   = score_gene('GATA6', preimpl_cells)
gata4   = score_gene('GATA4', preimpl_cells)
cdx2    = score_gene('CDX2',  preimpl_cells)
klf17   = score_gene('KLF17', preimpl_cells)

# Thresholds (log2 scale): use 75th percentile of positive cells
th_nanog = np.percentile(nanog[nanog > 0], 25) if (nanog > 0).sum() > 5 else 1.0
th_gata6 = np.percentile(gata6[gata6 > 0], 25) if (gata6 > 0).sum() > 5 else 1.0
th_gata4 = np.percentile(gata4[gata4 > 0], 25) if (gata4 > 0).sum() > 5 else 0.5
th_cdx2  = np.percentile(cdx2[cdx2 > 0],  25) if (cdx2 > 0).sum() > 5 else 0.5

print(f"Thresholds: NANOG={th_nanog:.2f}, GATA6={th_gata6:.2f}, GATA4={th_gata4:.2f}, CDX2={th_cdx2:.2f}")

nanog_pos = nanog >= th_nanog
gata6_pos = gata6 >= th_gata6
gata4_pos = gata4 >= th_gata4
cdx2_pos  = cdx2  >= th_cdx2

# Priority: TE > PE/Hypo > ICM (NANOG+GATA6+) > EPI (NANOG+GATA6-)
lineage = np.full(len(preimpl_cells), 'Unresolved', dtype=object)
lineage[nanog_pos & ~gata6_pos & ~cdx2_pos]               = 'EPI'
lineage[nanog_pos & gata6_pos & ~cdx2_pos]                = 'ICM'
lineage[gata6_pos & gata4_pos & ~nanog_pos & ~cdx2_pos]   = 'Hypo'
lineage[cdx2_pos]                                          = 'TE'

umap_df['lineage'] = lineage

from collections import Counter
lin_counts = Counter(lineage)
print("Lineage annotation:")
for k, n in sorted(lin_counts.items()):
    print(f"  {k}: {n} cells")

# ── 7. mTOR / quiescence signature ────────────────────────────────────────────
mtor_genes = ['FLCN', 'EIF4EBP1', 'MYC', 'MTOR', 'CDKN1A', 'TSC2']
ribo_genes  = [g for g in log_expr.index if
               (g.startswith('RPS') and not g.startswith('RPS6K')) or
               (g.startswith('RPL') and not g.startswith('RPLK'))]
print(f"Ribosomal genes available: {len(ribo_genes)}")

icm_mask   = lineage == 'ICM'
epi_mask   = lineage == 'EPI'
hypo_mask  = lineage == 'Hypo'
te_mask    = lineage == 'TE'

# Library size (raw RPKM sum per cell)
lib_size = expr.sum(axis=0).values  # per cell
lib_df = pd.DataFrame({'lib': lib_size, 'lineage': lineage}, index=preimpl_cells)

# Mean ribosomal gene expression per lineage
ribo_available = [g for g in ribo_genes if g in log_expr.index]
ribo_mean_icm  = log_expr.loc[ribo_available, np.array(preimpl_cells)[icm_mask]].values.mean(axis=1) if icm_mask.sum() else np.array([])
ribo_mean_epi  = log_expr.loc[ribo_available, np.array(preimpl_cells)[epi_mask]].values.mean(axis=1) if epi_mask.sum() else np.array([])

print(f"\n=== mTOR signature in cynomolgus ICM vs EPI ===")
print(f"ICM n={icm_mask.sum()}, EPI n={epi_mask.sum()}, Hypo n={hypo_mask.sum()}, TE n={te_mask.sum()}")

for gene in mtor_genes:
    if gene not in log_expr.index:
        print(f"  {gene}: NOT FOUND")
        continue
    vals = log_expr.loc[gene, preimpl_cells].values
    icm_vals = vals[icm_mask]
    epi_vals = vals[epi_mask]
    if icm_mask.sum() >= 3 and epi_mask.sum() >= 3:
        fc = np.median(icm_vals) - np.median(epi_vals)
        if icm_mask.sum() > 0 and epi_mask.sum() > 0:
            stat, p = mannwhitneyu(icm_vals, epi_vals, alternative='two-sided')
        else:
            p = np.nan
        print(f"  {gene}: ICM median={np.median(icm_vals):.2f}, EPI median={np.median(epi_vals):.2f}, log2FC={fc:+.2f}, p={p:.3g}")

# Ribosomal enrichment Fisher's exact test (ICM vs EPI)
if icm_mask.sum() >= 3 and epi_mask.sum() >= 3:
    all_gene_fc = []
    ribo_fc = []
    for gene in log_expr.index:
        icm_v = log_expr.loc[gene, np.array(preimpl_cells)[icm_mask]].values
        epi_v = log_expr.loc[gene, np.array(preimpl_cells)[epi_mask]].values
        fc = np.median(icm_v) - np.median(epi_v)
        all_gene_fc.append((gene, fc))

    all_fc_arr = np.array([x[1] for x in all_gene_fc])
    ribo_idx   = [i for i,(g,_) in enumerate(all_gene_fc) if g in ribo_available]
    ribo_fc_arr = all_fc_arr[ribo_idx]

    # Fraction significantly down (FC < -0.5) across all genes
    down_all  = (all_fc_arr < -0.5).sum()
    down_ribo = (ribo_fc_arr < -0.5).sum()
    n_all  = len(all_fc_arr)
    n_ribo = len(ribo_fc_arr)

    table = [[down_ribo, n_ribo - down_ribo],
             [down_all - down_ribo, n_all - n_ribo - (down_all - down_ribo)]]
    _, p_fisher = fisher_exact(table, alternative='greater')
    print(f"\nRibosomal enrichment (ICM<EPI, FC<-0.5): {down_ribo}/{n_ribo} ({100*down_ribo/n_ribo:.1f}%) vs background {down_all}/{n_all} ({100*down_all/n_all:.1f}%)")
    print(f"Fisher's exact p={p_fisher:.3g}")

# ── 8. Figures ────────────────────────────────────────────────────────────────
STAGE_COLORS = {'E06': '#4575b4', 'E07': '#74add1', 'E08': '#fdae61', 'E09': '#d73027'}
LIN_COLORS   = {'ICM': '#e41a1c', 'EPI': '#377eb8', 'Hypo': '#ff7f00', 'TE': '#4daf4a', 'Unresolved': '#999999'}

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Cynomolgus Macaque Preimplantation Embryo\n(Nakamura et al. 2017, GSE74767 — E06–E09, n=193)',
             fontsize=13, fontweight='bold')

u1, u2 = umap_df['UMAP1'].values, umap_df['UMAP2'].values

# Panel A: UMAP colored by stage
ax = axes[0, 0]
for stg, col in STAGE_COLORS.items():
    mask = stages.values == stg
    ax.scatter(u1[mask], u2[mask], c=col, s=30, alpha=0.8, label=stg, linewidths=0)
ax.set_title('A. Developmental stage', fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(fontsize=8, framealpha=0.7)
ax.set_aspect('equal', 'box')

# Panel B: UMAP colored by lineage
ax = axes[0, 1]
for lin, col in LIN_COLORS.items():
    mask = np.array(lineage) == lin
    if mask.sum() > 0:
        ax.scatter(u1[mask], u2[mask], c=col, s=35, alpha=0.85, label=f'{lin} (n={mask.sum()})', linewidths=0)
ax.set_title('B. Annotated lineage', fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(fontsize=8, framealpha=0.7)
ax.set_aspect('equal', 'box')

# Panel C: NANOG/GATA6 expression on UMAP
ax = axes[0, 2]
# Composite: colour by log2(NANOG) - log2(GATA6) to show ICM balance
nanog_vals = score_gene('NANOG', preimpl_cells)
gata6_vals = score_gene('GATA6', preimpl_cells)
score = nanog_vals - gata6_vals
sc = ax.scatter(u1, u2, c=score, cmap='RdBu_r', s=30, alpha=0.85, linewidths=0,
                vmin=-3, vmax=3)
plt.colorbar(sc, ax=ax, label='NANOG − GATA6 (log2)', shrink=0.8)
ax.set_title('C. NANOG vs GATA6 balance', fontweight='bold')
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.set_aspect('equal', 'box')

# Panel D: mTOR gene expression comparison (ICM vs EPI vs Hypo)
ax = axes[1, 0]
genes_to_plot = [g for g in ['FLCN','EIF4EBP1','MYC','CDKN1A','MTOR','TSC2'] if g in log_expr.index]
fc_icm_vs_epi  = []
fc_icm_vs_hypo = []
for gene in genes_to_plot:
    vals = log_expr.loc[gene, preimpl_cells].values
    icm_v  = vals[icm_mask]
    epi_v  = vals[epi_mask]
    hypo_v = vals[hypo_mask]
    fc_icm_vs_epi.append(np.median(icm_v) - np.median(epi_v) if epi_mask.sum() >= 3 else np.nan)
    fc_icm_vs_hypo.append(np.median(icm_v) - np.median(hypo_v) if hypo_mask.sum() >= 3 else np.nan)

x = np.arange(len(genes_to_plot))
w = 0.35
ax.bar(x - w/2, fc_icm_vs_epi,  width=w, color='#377eb8', label='ICM vs EPI',  alpha=0.85)
ax.bar(x + w/2, fc_icm_vs_hypo, width=w, color='#ff7f00', label='ICM vs Hypo', alpha=0.85)
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.set_xticks(x); ax.set_xticklabels(genes_to_plot, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('log₂FC (ICM vs other)')
ax.set_title('D. mTOR-quiescence genes: ICM vs EPI/Hypo', fontweight='bold')
ax.legend(fontsize=8)

# Panel E: Library size comparison
ax = axes[1, 1]
lin_order = ['ICM', 'EPI', 'Hypo', 'TE']
pos_data = []
for lin in lin_order:
    mask_lin = np.array(lineage) == lin
    if mask_lin.sum() > 0:
        pos_data.append(lib_size[mask_lin])
    else:
        pos_data.append(np.array([]))

parts = ax.violinplot([d for d in pos_data if len(d) > 0],
                      positions=[i for i, d in enumerate(pos_data) if len(d) > 0],
                      showmedians=True, showextrema=False)
valid_lins = [lin_order[i] for i, d in enumerate(pos_data) if len(d) > 0]
for i, (pc, lin) in enumerate(zip(parts['bodies'], valid_lins)):
    pc.set_facecolor(LIN_COLORS.get(lin, '#999'))
    pc.set_alpha(0.7)
ax.set_xticks(range(len(valid_lins)))
ax.set_xticklabels([f'{l}\n(n={sum(np.array(lineage)==l)})' for l in valid_lins], fontsize=8)
ax.set_ylabel('Total expression (RPKM sum)')
ax.set_title('E. Library size by lineage', fontweight='bold')

# Panel F: Ribosomal gene mean ICM vs EPI
ax = axes[1, 2]
ribo_plot = ribo_available[:20] if len(ribo_available) > 20 else ribo_available
icm_ribo  = [log_expr.loc[g, np.array(preimpl_cells)[icm_mask]].mean() for g in ribo_plot] if icm_mask.sum() else []
epi_ribo  = [log_expr.loc[g, np.array(preimpl_cells)[epi_mask]].mean() for g in ribo_plot] if epi_mask.sum() else []

if len(icm_ribo) > 0 and len(epi_ribo) > 0:
    ax.scatter(epi_ribo, icm_ribo, c='#333333', s=20, alpha=0.7)
    mn = min(min(icm_ribo), min(epi_ribo))
    mx = max(max(icm_ribo), max(epi_ribo))
    ax.plot([mn, mx], [mn, mx], 'r--', lw=1, label='y=x')
    ax.set_xlabel('Ribosomal gene mean (EPI)')
    ax.set_ylabel('Ribosomal gene mean (ICM)')
    frac_down = sum(np.array(icm_ribo) < np.array(epi_ribo)) / len(icm_ribo)
    ax.set_title(f'F. Ribosomal gene expression ICM vs EPI\n({frac_down*100:.0f}% lower in ICM, n={len(ribo_plot)} genes)', fontweight='bold')
    ax.legend(fontsize=8)

plt.tight_layout()
outpath = '/Users/adelosan/cESFW_analysis/plots/FigNHP_Cynomolgus_Preimplantation.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
print(f"\nSaved: {outpath}")
plt.close()

# ── 9. Print summary statistics ────────────────────────────────────────────────
print("\n=== CYNOMOLGUS SUMMARY ===")
print(f"Total preimplantation cells: {len(preimpl_cells)}")
print(f"Lineage annotation (marker-based):")
for lin in ['ICM', 'EPI', 'Hypo', 'TE', 'Unresolved']:
    n = (np.array(lineage) == lin).sum()
    print(f"  {lin}: {n} cells ({100*n/len(lineage):.1f}%)")

# Library size ratio ICM vs EPI
if icm_mask.sum() >= 3 and epi_mask.sum() >= 3:
    lib_icm = lib_size[icm_mask]
    lib_epi = lib_size[epi_mask]
    _, p_lib = mannwhitneyu(lib_icm, lib_epi, alternative='two-sided')
    ratio = np.median(lib_icm) / np.median(lib_epi)
    print(f"\nLibrary size: ICM/EPI ratio = {ratio:.3f}, Mann-Whitney p = {p_lib:.3g}")

# Ribosomal gene depletion
if icm_mask.sum() >= 3 and epi_mask.sum() >= 3:
    n_down_ribo = sum(
        log_expr.loc[g, np.array(preimpl_cells)[icm_mask]].median() <
        log_expr.loc[g, np.array(preimpl_cells)[epi_mask]].median()
        for g in ribo_available
    )
    print(f"Ribosomal genes lower in ICM vs EPI: {n_down_ribo}/{len(ribo_available)} ({100*n_down_ribo/len(ribo_available):.1f}%)")
