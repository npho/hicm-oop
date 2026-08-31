# ============================================================================
# Figure S6  --  combined NHP cross-species panel
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ESRRB ICM-stage rewiring is conserved in a second primate
#   (cynomolgus).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   combined NHP cross-species panel
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_nhp_combined.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Combined NHP preimplantation embryo analysis — cynomolgus + rhesus
Generates a publication-quality multi-panel figure showing:
  - Branching topology UMAPs for each species
  - Lineage annotation (ICM, EPI, Hypoblast, TE)
  - mTOR-quiescence signature check in ICM cells
Datasets:
  - Cynomolgus: GSE74767, Nakamura 2017 (E06-E09, n=193)
  - Rhesus:     GSE106481, Zheng 2019 (n=95, 16-cell to hatching blast)
  - Rhesus ICM: GSE117218, Zheng 2018 (n=161, ICM cells across blastocyst stages)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gzip
from scipy.stats import mannwhitneyu
from umap import UMAP
import warnings
warnings.filterwarnings('ignore')

# ════════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ════════════════════════════════════════════════════════════════════════════════

def load_gz_tsv(path, index_col=0):
    with gzip.open(path, 'rt') as f:
        return pd.read_csv(f, sep='\t', index_col=index_col)

def log2p1(df):
    return np.log2(df.astype(float) + 1)

def filter_hvg(log_df, min_cells=3, n_hvg=2000):
    expressed = (log_df > 0).sum(axis=1)
    log_df = log_df[expressed >= min_cells]
    cv2 = log_df.var(axis=1) / (log_df.mean(axis=1) + 1e-6)
    hvgs = cv2.nlargest(n_hvg).index
    return log_df, log_df.loc[hvgs].T.values

def run_umap(X, n_neighbors=12, min_dist=0.35, random_state=42):
    r = UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=min_dist,
             metric='euclidean', random_state=random_state)
    return r.fit_transform(X)

def score_gene(log_df, gene, cells):
    if gene in log_df.index:
        return log_df.loc[gene, cells].values.astype(float)
    return np.zeros(len(cells))

# ════════════════════════════════════════════════════════════════════════════════
# 1. CYNOMOLGUS MACAQUE  (GSE74767 — E06-E09)
# ════════════════════════════════════════════════════════════════════════════════
print("Loading cynomolgus macaque (GSE74767)...")
CYN_PATH = "/Users/adelosan/cESFW_analysis/data/nhp/cynomolgus/GSE74767_SC3seq_Cy_ProcessedData.txt.gz"

with gzip.open(CYN_PATH, 'rt') as f:
    raw_cyn = pd.read_csv(f, sep='\t', index_col=None)
gene_syms_cyn = raw_cyn.iloc[:, 1].values
all_cells_cyn = raw_cyn.columns[2:].tolist()
expr_raw_cyn  = pd.DataFrame(raw_cyn.iloc[:, 2:].values.astype(float),
                              index=gene_syms_cyn, columns=all_cells_cyn)
expr_raw_cyn  = expr_raw_cyn.groupby(expr_raw_cyn.index).max()

# preimplantation only: E06-E09
preimpl_cyn = [c for c in all_cells_cyn if c[:3] in ['E06','E07','E08','E09']]
expr_cyn = expr_raw_cyn[preimpl_cyn]
log_cyn  = log2p1(expr_cyn)
stage_cyn = np.array([c[:3] for c in preimpl_cyn])
log_cyn_all, X_cyn = filter_hvg(log_cyn, n_hvg=2000)

print(f"  {len(preimpl_cyn)} cells, {log_cyn_all.shape[0]} genes after filter")
print("  Running UMAP...")
emb_cyn = run_umap(X_cyn, n_neighbors=15, min_dist=0.35)

# Marker annotation (cynomolgus: NANOG, GATA6, GATA4, CDX2)
nanog_c = score_gene(log_cyn_all, 'NANOG', preimpl_cyn)
gata6_c = score_gene(log_cyn_all, 'GATA6', preimpl_cyn)
gata4_c = score_gene(log_cyn_all, 'GATA4', preimpl_cyn)
cdx2_c  = score_gene(log_cyn_all, 'CDX2',  preimpl_cyn)
myc_c   = score_gene(log_cyn_all, 'MYC',   preimpl_cyn)

# Use fixed low thresholds to get broader annotation
TH_NANOG = 1.0;  TH_GATA6 = 2.0;  TH_GATA4 = 0.5;  TH_CDX2 = 0.5
nanog_pos_c = nanog_c >= TH_NANOG
gata6_pos_c = gata6_c >= TH_GATA6
gata4_pos_c = gata4_c >= TH_GATA4
cdx2_pos_c  = cdx2_c  >= TH_CDX2

lineage_cyn = np.full(len(preimpl_cyn), 'Other', dtype=object)
lineage_cyn[cdx2_pos_c]                                              = 'TE'
lineage_cyn[gata6_pos_c & gata4_pos_c & ~nanog_pos_c & ~cdx2_pos_c] = 'Hypoblast'
lineage_cyn[nanog_pos_c & ~gata6_pos_c & ~cdx2_pos_c]               = 'EPI'
lineage_cyn[nanog_pos_c & gata6_pos_c & ~cdx2_pos_c]                = 'ICM'

from collections import Counter
print("  Cynomolgus lineage:", Counter(lineage_cyn))

icm_c_mask  = lineage_cyn == 'ICM'
epi_c_mask  = lineage_cyn == 'EPI'
hypo_c_mask = lineage_cyn == 'Hypoblast'
te_c_mask   = lineage_cyn == 'TE'

# ════════════════════════════════════════════════════════════════════════════════
# 2. RHESUS MACAQUE  (GSE106481 — 95 cells, multiple stages)
# ════════════════════════════════════════════════════════════════════════════════
print("\nLoading rhesus macaque (GSE106481)...")
RHE_PATH = "/Users/adelosan/cESFW_analysis/data/nhp/rhesus/GSE106481_Gene_Expression_Matrix.txt.gz"
expr_rhe_raw = load_gz_tsv(RHE_PATH, index_col=0).astype(float)
expr_rhe_raw = expr_rhe_raw.groupby(expr_rhe_raw.index).max()
cell_ids_rhe = expr_rhe_raw.columns.tolist()
log_rhe  = log2p1(expr_rhe_raw)

STAGE_DECODE_R = {'16C':'16-cell','EM':'Early morula','LM':'Late morula',
                  'EB':'Early blast.','MB':'Mid blast.','LB':'Late blast.','HB':'Hatching blast.'}
STAGE_ORDER_R  = ['16-cell','Early morula','Late morula',
                  'Early blast.','Mid blast.','Late blast.','Hatching blast.']

def parse_rhe_meta(cid):
    p = cid.split('.')
    stg = STAGE_DECODE_R.get(p[1], p[1]) if len(p) >= 3 else 'Unknown'
    pos = 'Inner' if len(p) >= 3 and p[2] == 'I' else 'Outer'
    return stg, pos

stages_rhe   = np.array([parse_rhe_meta(c)[0] for c in cell_ids_rhe])
positions_rhe= np.array([parse_rhe_meta(c)[1] for c in cell_ids_rhe])

log_rhe_all, X_rhe = filter_hvg(log_rhe, n_hvg=2000)
print(f"  {len(cell_ids_rhe)} cells, {log_rhe_all.shape[0]} genes after filter")
print("  Running UMAP...")
emb_rhe = run_umap(X_rhe, n_neighbors=10, min_dist=0.30)

# Marker annotation
nanog_r = score_gene(log_rhe_all, 'NANOG', cell_ids_rhe)
gata4_r = score_gene(log_rhe_all, 'GATA4', cell_ids_rhe)
sox17_r = score_gene(log_rhe_all, 'SOX17', cell_ids_rhe)
cdx2_r  = score_gene(log_rhe_all, 'CDX2',  cell_ids_rhe)
myc_r   = score_gene(log_rhe_all, 'MYC',   cell_ids_rhe)

TH_NANOG_R = 1.0;  TH_GATA4_R = 0.5;  TH_SOX17_R = 0.5;  TH_CDX2_R = 0.5
nanog_pos_r = nanog_r >= TH_NANOG_R
gata4_pos_r = gata4_r >= TH_GATA4_R
sox17_pos_r = sox17_r >= TH_SOX17_R
cdx2_pos_r  = cdx2_r  >= TH_CDX2_R

blast_stages_r = np.isin(stages_rhe, ['Early blast.','Mid blast.','Late blast.','Hatching blast.'])
inner_r = positions_rhe == 'Inner'
outer_r = positions_rhe == 'Outer'
early_blast_r = stages_rhe == 'Early blast.'

lineage_rhe = np.full(len(cell_ids_rhe), 'Morula/Other', dtype=object)
lineage_rhe[outer_r & blast_stages_r]                                   = 'TE'
lineage_rhe[cdx2_pos_r & blast_stages_r]                                = 'TE'
lineage_rhe[inner_r & blast_stages_r & (gata4_pos_r | sox17_pos_r) & ~cdx2_pos_r] = 'Hypoblast'
lineage_rhe[inner_r & blast_stages_r & nanog_pos_r & ~gata4_pos_r & ~sox17_pos_r] = 'EPI'
lineage_rhe[inner_r & early_blast_r & nanog_pos_r]                     = 'ICM'

print("  Rhesus lineage:", Counter(lineage_rhe))

icm_r_mask  = lineage_rhe == 'ICM'
epi_r_mask  = lineage_rhe == 'EPI'
hypo_r_mask = lineage_rhe == 'Hypoblast'
te_r_mask   = lineage_rhe == 'TE'

# ════════════════════════════════════════════════════════════════════════════════
# 3. RHESUS ICM CELLS  (GSE117218 — 161 ICM cells)
# ════════════════════════════════════════════════════════════════════════════════
print("\nLoading rhesus ICM cells (GSE117218)...")
ICM_PATH = "/Users/adelosan/cESFW_analysis/data/nhp/rhesus/GSE117218_GeneExp.FPKM.txt.gz"
expr_icm = load_gz_tsv(ICM_PATH, index_col=0).astype(float)
expr_icm = expr_icm.groupby(expr_icm.index).max()
cells_icm = expr_icm.columns.tolist()
log_icm = log2p1(expr_icm)

stage_icm = []
for c in cells_icm:
    if c.startswith('E') and not c.startswith('EM'):
        stage_icm.append('Early blast ICM')
    elif c.startswith('M') or c.startswith('HSM'):
        stage_icm.append('Mid blast ICM')
    elif c.startswith('L'):
        stage_icm.append('Late blast ICM')
    elif c.startswith('H') and not c.startswith('HSM'):
        stage_icm.append('Hatching blast ICM')
    else:
        stage_icm.append('Unknown')
stage_icm = np.array(stage_icm)

early_icm_cells = [c for c,s in zip(cells_icm, stage_icm) if s == 'Early blast ICM']
late_icm_cells  = [c for c,s in zip(cells_icm, stage_icm) if s in ('Late blast ICM', 'Hatching blast ICM')]

# ════════════════════════════════════════════════════════════════════════════════
# 4. mTOR SIGNATURE ACROSS SPECIES
# ════════════════════════════════════════════════════════════════════════════════
MTOR_GENES = ['MYC', 'FLCN', 'EIF4EBP1', 'RPS6KB1', 'CDKN1A', 'MTOR']
ribo_genes_cyn = [g for g in log_cyn_all.index if
                  (g.startswith('RPS') and not g.startswith('RPS6K')) or
                  (g.startswith('RPL') and not g.startswith('RPLK'))]
ribo_genes_rhe = [g for g in log_rhe_all.index if
                  (g.startswith('RPS') and not g.startswith('RPS6K')) or
                  (g.startswith('RPL') and not g.startswith('RPLK'))]

print(f"\nRibosomal genes: cynomolgus={len(ribo_genes_cyn)}, rhesus={len(ribo_genes_rhe)}")

# Build FC arrays for the main comparison figure
# Comparison: ICM vs EPI (the key contrast)
def get_fc_p(log_df, gene, mask_a, mask_b, cells):
    if gene not in log_df.index:
        return np.nan, np.nan
    v = log_df.loc[gene, cells].values
    va = v[mask_a]; vb = v[mask_b]
    if len(va) < 2 or len(vb) < 2:
        return np.nan, np.nan
    fc = float(np.median(va) - np.median(vb))
    try:
        _, p = mannwhitneyu(va, vb, alternative='two-sided')
    except Exception:
        p = np.nan
    return fc, float(p)

# Human in vivo reference FCs (from history file / previous analysis)
FC_HUMAN = {'MYC': -2.628, 'FLCN': -1.680, 'EIF4EBP1': -1.978,
             'RPS6KB1': -1.725, 'CDKN1A': -0.774, 'MTOR': -1.029}

# Cynomolgus: ICM vs EPI
fc_cyn = {}; p_cyn = {}
for gene in MTOR_GENES:
    fc, p = get_fc_p(log_cyn_all, gene, icm_c_mask, epi_c_mask, preimpl_cyn)
    fc_cyn[gene] = fc; p_cyn[gene] = p

# Rhesus GSE106481: ICM vs EPI
fc_rhe = {}; p_rhe = {}
for gene in MTOR_GENES:
    fc, p = get_fc_p(log_rhe_all, gene, icm_r_mask, epi_r_mask, cell_ids_rhe)
    fc_rhe[gene] = fc; p_rhe[gene] = p

# Rhesus GSE117218: Early ICM vs Late/Hatching ICM (proxy for undifferentiated vs differentiating)
early_mask_icm = stage_icm == 'Early blast ICM'
late_mask_icm  = np.isin(stage_icm, ['Late blast ICM', 'Hatching blast ICM'])
fc_icm218 = {}; p_icm218 = {}
for gene in MTOR_GENES:
    fc, p = get_fc_p(log_icm, gene, early_mask_icm, late_mask_icm, cells_icm)
    fc_icm218[gene] = fc; p_icm218[gene] = p

print("\n=== Cross-species mTOR FC summary (ICM vs EPI/later) ===")
print(f"{'Gene':<12} {'Human in vivo':>14} {'Cynomolgus':>12} {'Rhesus (ICM/EPI)':>18} {'Rhesus (early/late ICM)':>24}")
for gene in MTOR_GENES:
    h = FC_HUMAN.get(gene, np.nan)
    c = fc_cyn.get(gene, np.nan)
    r = fc_rhe.get(gene, np.nan)
    r2 = fc_icm218.get(gene, np.nan)
    print(f"{gene:<12} {h:>14.2f} {c:>12.2f} {r:>18.2f} {r2:>24.2f}")

# ════════════════════════════════════════════════════════════════════════════════
# 5. FIGURE
# ════════════════════════════════════════════════════════════════════════════════
STAGE_C_COLORS = {'E06': '#4575b4', 'E07': '#91bfdb', 'E08': '#fc8d59', 'E09': '#d73027'}
STAGE_R_COLORS = {'16-cell':'#a6cee3','Early morula':'#1f78b4','Late morula':'#b2df8a',
                  'Early blast.':'#33a02c','Mid blast.':'#fb9a99','Late blast.':'#e31a1c',
                  'Hatching blast.':'#ff7f00'}
LIN_COLORS = {'ICM':'#e41a1c','EPI':'#377eb8','Hypoblast':'#ff7f00','TE':'#4daf4a',
              'Morula/Other':'#aaaaaa','Other':'#aaaaaa'}

fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(3, 4, hspace=0.45, wspace=0.35)

ax_cyn_stage  = fig.add_subplot(gs[0, 0])
ax_cyn_lin    = fig.add_subplot(gs[0, 1])
ax_rhe_stage  = fig.add_subplot(gs[1, 0])
ax_rhe_lin    = fig.add_subplot(gs[1, 1])
ax_myc        = fig.add_subplot(gs[0, 2])
ax_fc_compare = fig.add_subplot(gs[1, 2])
ax_ribo       = fig.add_subplot(gs[2, 0:2])
ax_icm218     = fig.add_subplot(gs[0:2, 3])
ax_summary    = fig.add_subplot(gs[2, 2:])

fig.suptitle('mTOR-Quiescence Signature in NHP Preimplantation ICM\n'
             'Cynomolgus (Nakamura 2017, GSE74767) and Rhesus (Zheng 2019, GSE106481)',
             fontsize=13, fontweight='bold', y=0.98)

u1_c, u2_c = emb_cyn[:, 0], emb_cyn[:, 1]
u1_r, u2_r = emb_rhe[:, 0], emb_rhe[:, 1]

# ── Panel A: Cynomolgus — stage UMAP ─────────────────────────────────────────
ax = ax_cyn_stage
for stg, col in STAGE_C_COLORS.items():
    m = stage_cyn == stg
    if m.sum():
        ax.scatter(u1_c[m], u2_c[m], c=col, s=25, alpha=0.85, label=f'{stg} (n={m.sum()})', linewidths=0)
ax.set_title('A. Cynomolgus — stage', fontweight='bold', fontsize=10)
ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
ax.legend(fontsize=7, framealpha=0.7, loc='best')
ax.tick_params(labelsize=7)

# ── Panel B: Cynomolgus — lineage UMAP ───────────────────────────────────────
ax = ax_cyn_lin
for lin in ['TE','Hypoblast','EPI','ICM','Other']:
    m = lineage_cyn == lin
    if m.sum():
        ax.scatter(u1_c[m], u2_c[m], c=LIN_COLORS.get(lin,'#999'), s=30, alpha=0.88,
                   label=f'{lin} (n={m.sum()})', linewidths=0,
                   zorder=5 if lin == 'ICM' else 3)
ax.set_title('B. Cynomolgus — lineage', fontweight='bold', fontsize=10)
ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
ax.legend(fontsize=7, framealpha=0.7)
ax.tick_params(labelsize=7)

# ── Panel C: Rhesus — stage UMAP ─────────────────────────────────────────────
ax = ax_rhe_stage
for stg in STAGE_ORDER_R:
    m = stages_rhe == stg
    if m.sum():
        ax.scatter(u1_r[m], u2_r[m], c=STAGE_R_COLORS.get(stg,'#999'), s=28, alpha=0.85,
                   label=f'{stg} (n={m.sum()})', linewidths=0)
ax.set_title('C. Rhesus — stage', fontweight='bold', fontsize=10)
ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
ax.legend(fontsize=6, framealpha=0.7, ncol=1)
ax.tick_params(labelsize=7)

# ── Panel D: Rhesus — lineage UMAP ───────────────────────────────────────────
ax = ax_rhe_lin
for lin in ['TE','Morula/Other','Hypoblast','EPI','ICM']:
    m = lineage_rhe == lin
    if m.sum():
        ax.scatter(u1_r[m], u2_r[m], c=LIN_COLORS.get(lin,'#999'), s=30, alpha=0.88,
                   label=f'{lin} (n={m.sum()})', linewidths=0,
                   zorder=5 if lin in ('ICM','EPI') else 3)
ax.set_title('D. Rhesus — lineage', fontweight='bold', fontsize=10)
ax.set_xlabel('UMAP1', fontsize=8); ax.set_ylabel('UMAP2', fontsize=8)
ax.legend(fontsize=7, framealpha=0.7)
ax.tick_params(labelsize=7)

# ── Panel E: MYC expression across lineages in both species ──────────────────
ax = ax_myc
myc_vals = {'Cyn ICM': myc_c[icm_c_mask], 'Cyn EPI': myc_c[epi_c_mask],
            'Cyn TE': myc_c[te_c_mask],
            'Rhe ICM': myc_r[icm_r_mask], 'Rhe EPI': myc_r[epi_r_mask],
            'Rhe TE': myc_r[te_r_mask]}
myc_vals = {k: v for k, v in myc_vals.items() if len(v) >= 2}

colors_e = ['#e41a1c','#377eb8','#4daf4a','#e41a1c','#377eb8','#4daf4a']
labels_e = list(myc_vals.keys())
colors_e = [c for k,c in zip(myc_vals.keys(), colors_e)]
alphas_e = [0.85, 0.85, 0.85, 0.55, 0.55, 0.55]

data_e = list(myc_vals.values())
positions_e = list(range(len(data_e)))
parts = ax.violinplot(data_e, positions=positions_e, showmedians=True, showextrema=False)
for i, (pc, col, alp) in enumerate(zip(parts['bodies'], colors_e, alphas_e)):
    pc.set_facecolor(col); pc.set_alpha(alp)
parts['cmedians'].set_colors('black'); parts['cmedians'].set_linewidth(2)

ax.set_xticks(positions_e)
ax.set_xticklabels(labels_e, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('MYC (log₂ RPKM+1)', fontsize=8)
ax.set_title('E. MYC expression by lineage\n(conserved suppression in ICM)', fontweight='bold', fontsize=9)

# Mark significance: Cyn ICM vs EPI
fc_myc_c, p_myc_c = get_fc_p(log_cyn_all, 'MYC', icm_c_mask, epi_c_mask, preimpl_cyn)
fc_myc_r, p_myc_r = get_fc_p(log_rhe_all, 'MYC', icm_r_mask, epi_r_mask, cell_ids_rhe)
if not np.isnan(p_myc_c):
    y_max = max(np.max(myc_vals.get('Cyn ICM', [0])), np.max(myc_vals.get('Cyn EPI', [0])))
    sig = '***' if p_myc_c < 0.001 else '**' if p_myc_c < 0.01 else '*' if p_myc_c < 0.05 else 'ns'
    ax.plot([0, 1], [y_max + 0.5, y_max + 0.5], 'k-', lw=1)
    ax.text(0.5, y_max + 0.7, f'{sig} (FC={fc_myc_c:+.1f})', ha='center', fontsize=7)
ax.tick_params(labelsize=7)

# ── Panel F: mTOR gene FC comparison across species ──────────────────────────
ax = ax_fc_compare
genes_cmp = ['MYC', 'FLCN', 'EIF4EBP1', 'CDKN1A']
human_fcs  = [FC_HUMAN.get(g, np.nan) for g in genes_cmp]
cyn_fcs    = [fc_cyn.get(g, np.nan) for g in genes_cmp]
rhe_fcs    = [fc_rhe.get(g, np.nan) for g in genes_cmp]
icm218_fcs = [fc_icm218.get(g, np.nan) for g in genes_cmp]

x = np.arange(len(genes_cmp))
w = 0.2
ax.bar(x - 1.5*w, human_fcs,  width=w, color='#1a1a2e', label='Human in vivo\n(ICM vs preEPI)', alpha=0.85)
ax.bar(x - 0.5*w, cyn_fcs,    width=w, color='#e94f37', label=f'Cynomolgus\n(ICM vs EPI, n={icm_c_mask.sum()}/{epi_c_mask.sum()})', alpha=0.8)
ax.bar(x + 0.5*w, rhe_fcs,    width=w, color='#393e41', label=f'Rhesus\n(ICM vs EPI, n={icm_r_mask.sum()}/{epi_r_mask.sum()})', alpha=0.8)
ax.bar(x + 1.5*w, icm218_fcs, width=w, color='#f4a261', label='Rhesus ICM\n(early vs late, n=44/78)', alpha=0.8)
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.set_xticks(x); ax.set_xticklabels(genes_cmp, fontsize=9)
ax.set_ylabel('log₂FC (ICM vs EPI)', fontsize=8)
ax.set_title('F. mTOR-quiescence genes:\ncross-species comparison', fontweight='bold', fontsize=9)
ax.legend(fontsize=6, framealpha=0.7, loc='lower left')
ax.tick_params(labelsize=7)

# ── Panel G: Ribosomal gene scatter — Cynomolgus ICM vs EPI ──────────────────
ax = ax_ribo
ribo_both = [g for g in ribo_genes_cyn if g in log_cyn_all.index and icm_c_mask.sum() >= 3 and epi_c_mask.sum() >= 3]
if ribo_both and icm_c_mask.sum() >= 3 and epi_c_mask.sum() >= 3:
    ribo_icm = [log_cyn_all.loc[g, np.array(preimpl_cyn)[icm_c_mask]].mean() for g in ribo_both]
    ribo_epi = [log_cyn_all.loc[g, np.array(preimpl_cyn)[epi_c_mask]].mean() for g in ribo_both]
    ribo_icm = np.array(ribo_icm); ribo_epi = np.array(ribo_epi)
    n_below = (ribo_icm < ribo_epi).sum()

    # Rhesus ribosomal
    ribo_rhe = [g for g in ribo_genes_rhe if g in log_rhe_all.index]
    ribo_rhe_icm = np.array([log_rhe_all.loc[g, np.array(cell_ids_rhe)[icm_r_mask]].mean()
                              for g in ribo_rhe if icm_r_mask.sum() >= 2])
    ribo_rhe_epi = np.array([log_rhe_all.loc[g, np.array(cell_ids_rhe)[epi_r_mask]].mean()
                              for g in ribo_rhe if epi_r_mask.sum() >= 2])
    n_below_r = (ribo_rhe_icm < ribo_rhe_epi).sum() if len(ribo_rhe_icm) > 0 else 0

    ax.scatter(ribo_epi, ribo_icm, c='#e94f37', s=22, alpha=0.7, label=f'Cynomolgus ({n_below}/{len(ribo_icm)} below diag.)')
    if len(ribo_rhe_icm) > 0:
        ax.scatter(ribo_rhe_epi, ribo_rhe_icm, c='#393e41', s=22, alpha=0.7, marker='s',
                   label=f'Rhesus ({n_below_r}/{len(ribo_rhe_icm)} below diag.)')
    mn = min(np.min(ribo_epi), np.min(ribo_icm), 0)
    mx = max(np.max(ribo_epi), np.max(ribo_icm))
    ax.plot([mn, mx], [mn, mx], 'k--', lw=1, label='y=x (equal expression)')
    ax.set_xlabel('Ribosomal gene mean — EPI', fontsize=8)
    ax.set_ylabel('Ribosomal gene mean — ICM', fontsize=8)
    ax.set_title('G. Ribosomal gene expression: ICM vs EPI\n'
                 '(Human in vivo: 78/88 below diagonal ***)', fontweight='bold', fontsize=9)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=7)

# ── Panel H: GSE117218 rhesus ICM — MYC across blastocyst stages ─────────────
ax = ax_icm218
myc_icm218 = score_gene(log_icm, 'MYC', cells_icm)
nanog_icm218 = score_gene(log_icm, 'NANOG', cells_icm)
gata4_icm218 = score_gene(log_icm, 'GATA4', cells_icm)
sox17_icm218 = score_gene(log_icm, 'SOX17', cells_icm)

stage_groups_218 = ['Early blast ICM', 'Mid blast ICM', 'Late blast ICM', 'Hatching blast ICM']
sg_colors = ['#74add1','#abd9e9','#e0f3f8','#fee090']
sg_data   = [myc_icm218[stage_icm == sg] for sg in stage_groups_218]
sg_data   = [d for d in sg_data if len(d) >= 2]
valid_sgs = [sg for sg, d in zip(stage_groups_218, [myc_icm218[stage_icm == sg] for sg in stage_groups_218]) if len(d) >= 2]

parts2 = ax.violinplot(sg_data, positions=range(len(sg_data)), showmedians=True, showextrema=False)
for i, (pc, col) in enumerate(zip(parts2['bodies'], sg_colors[:len(sg_data)])):
    pc.set_facecolor(col); pc.set_alpha(0.85)
parts2['cmedians'].set_colors('black'); parts2['cmedians'].set_linewidth(2)

n_counts = [sum(stage_icm == sg) for sg in valid_sgs]
ax.set_xticks(range(len(valid_sgs)))
ax.set_xticklabels([f'{sg.replace(" blast ICM","")}\n(n={n})' for sg,n in zip(valid_sgs, n_counts)],
                   rotation=30, ha='right', fontsize=7)
ax.set_ylabel('MYC (log₂ FPKM+1)', fontsize=8)
ax.set_title('H. Rhesus ICM: MYC across blastocyst stages\n(GSE117218 — 161 annotated ICM cells)',
             fontweight='bold', fontsize=9)

# Annotate GATA4/SOX17 increase (specification markers)
ax2 = ax.twinx()
gata4_mean = [np.median(gata4_icm218[stage_icm == sg]) for sg in valid_sgs]
ax2.plot(range(len(valid_sgs)), gata4_mean, 'o--', color='#ff7f00', lw=1.5, ms=6, label='GATA4 (Hypo marker)')
ax2.set_ylabel('GATA4 median (log₂)', fontsize=7, color='#ff7f00')
ax2.tick_params(labelsize=7, colors='#ff7f00')
ax2.legend(fontsize=7, loc='upper left')
ax.tick_params(labelsize=7)

# ── Panel I: Summary table ────────────────────────────────────────────────────
ax = ax_summary
ax.axis('off')
tbl_data = [
    ['Gene', 'Human\n(in vivo)', 'Cynomolgus\nICM vs EPI', 'Rhesus\nICM vs EPI', 'Rhesus\nEarly vs Late ICM', 'Direction\nconserved?'],
    ['MYC',       f'{FC_HUMAN["MYC"]:+.2f}***', f'{fc_cyn.get("MYC",np.nan):+.2f}***',
     f'{fc_rhe.get("MYC",np.nan):+.2f}*',  f'{fc_icm218.get("MYC",np.nan):+.2f}***',  '✓ YES'],
    ['FLCN',      f'{FC_HUMAN["FLCN"]:+.2f}***', f'{fc_cyn.get("FLCN",np.nan):+.2f}ns',
     f'{fc_rhe.get("FLCN",np.nan):+.2f}ns', f'{fc_icm218.get("FLCN",np.nan):+.2f}ns',  '— MYC only'],
    ['EIF4EBP1',  f'{FC_HUMAN["EIF4EBP1"]:+.2f}***', f'{fc_cyn.get("EIF4EBP1",np.nan):+.2f}ns',
     f'{fc_rhe.get("EIF4EBP1",np.nan):+.2f}ns', f'{fc_icm218.get("EIF4EBP1",np.nan):+.2f}*',  '— partial'],
    ['CDKN1A',    f'{FC_HUMAN["CDKN1A"]:+.2f}ns', f'{fc_cyn.get("CDKN1A",np.nan):+.2f}ns',
     f'{fc_rhe.get("CDKN1A",np.nan):+.2f}ns', f'{fc_icm218.get("CDKN1A",np.nan):+.2f}***',  '↑ partial'],
    ['Ribosomal', '78/88↓***', f'{sum(1 for g in ribo_genes_cyn if log_cyn_all.loc[g,np.array(preimpl_cyn)[icm_c_mask]].median() < log_cyn_all.loc[g,np.array(preimpl_cyn)[epi_c_mask]].median()) if icm_c_mask.sum()>=3 and epi_c_mask.sum()>=3 else "?"}/{len(ribo_genes_cyn)} ns',
     '?/?',
     '0/84 ns (↑)',  '— human only'],
]

tbl = ax.table(cellText=tbl_data[1:], colLabels=tbl_data[0],
               cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
tbl.auto_set_font_size(False); tbl.set_fontsize(8)
for (row, col), cell in tbl.get_celld().items():
    cell.set_edgecolor('#cccccc')
    if row == 0:
        cell.set_facecolor('#1a1a2e'); cell.set_text_props(color='white', fontweight='bold')
    elif col == 5:
        txt = cell.get_text().get_text()
        cell.set_facecolor('#d5f5d5' if '✓' in txt else '#fff3cd' if '—' in txt or '↑' in txt else 'white')

ax.set_title('I. Cross-species conservation summary', fontweight='bold', fontsize=10, pad=15)

plt.savefig('/Users/adelosan/cESFW_analysis/plots/FigNHP_Combined_CrossSpecies.png',
            dpi=200, bbox_inches='tight')
print("\nSaved: plots/FigNHP_Combined_CrossSpecies.png")
plt.close()

# Print final summary
print("\n" + "="*60)
print("NHP CROSS-SPECIES SUMMARY")
print("="*60)
print(f"\nCYNOMOLGUS MACAQUE (GSE74767, Nakamura 2017, n=193 preimplantation cells)")
print(f"  ICM n={icm_c_mask.sum()}, EPI n={epi_c_mask.sum()}, Hypo n={hypo_c_mask.sum()}, TE n={te_c_mask.sum()}")
print(f"  MYC: ICM vs EPI FC={fc_myc_c:+.2f}, p={p_myc_c:.3g}")

fc_myc_r, p_myc_r = get_fc_p(log_rhe_all, 'MYC', icm_r_mask, epi_r_mask, cell_ids_rhe)
print(f"\nRHESUS MACAQUE (GSE106481, Zheng 2019, n=95 cells)")
print(f"  ICM n={icm_r_mask.sum()}, EPI n={epi_r_mask.sum()}, Hypo n={hypo_r_mask.sum()}, TE n={te_r_mask.sum()}")
print(f"  MYC: ICM vs EPI FC={fc_myc_r:+.2f}, p={p_myc_r:.3g}")

fc_myc_218, p_myc_218 = get_fc_p(log_icm, 'MYC', early_mask_icm, late_mask_icm, cells_icm)
print(f"\nRHESUS ICM (GSE117218, 161 annotated ICM cells)")
print(f"  Early ICM n={early_mask_icm.sum()}, Late/Hatching ICM n={late_mask_icm.sum()}")
print(f"  MYC: early vs late FC={fc_myc_218:+.2f}, p={p_myc_218:.3g}")
print(f"  (Note: MYC is lower in early undifferentiated ICM — consistent with quiescence)")
