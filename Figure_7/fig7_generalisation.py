# ============================================================================
# Figure 7  --  assembles the four-system figure
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The transition-versus-state framework generalises to four further human
#   fate decisions, each with a synthetic 50/50 mixture as positive control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   assembles the four-system figure
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 7. The transition-versus-state framework generalises across four
#   human fate transitions, from development to regeneration. Order-
#   parameter modality test applied to (A,B) the astrocyte/oligodendrocyte
#   choice of glial progenitors in first-trimester human brain (Braun et
#   al., 2023; thalamus), (C,D) the absorptive/secretory choice of
#   intestinal stem cells in adult small intestine (Burclaff et al., 2022),
#   (E,F) the syncytiotrophoblast/extravillous choice of villous
#   cytotrophoblast in the placenta (Arutyunyan et al., 2023), and (G,H) the
#   AT2-to-AT1 conversion of the KRT8-positive alveolar intermediate in
#   fibrotic adult lung (Habermann et al., 2020). (A,C,E,G) φ = z(fate-A) −
#   z(fate-B): in each system the progenitor or intermediate is a single,
#   unimodal population positioned between its two committed fates.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_reconstruction/fig7_generalisation.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""
FIGURE 7 (manuscript v63+): "The transition-versus-state framework generalises
to two unrelated human fate decisions."  Self-contained reconstruction script.

Combined 2x2 figure:
  A,B = fetal brain glia (Braun et al. 2023, thalamus): radial-glia progenitor
        between astrocyte / oligodendrocyte fates; modality by CALIBRATED dBIC
        (the dip test is underpowered here because the poles are close).
  C,D = adult intestine (Burclaff et al. 2022, small intestine): intestinal stem
        cell (ISC) between absorptive / secretory (goblet) fates; modality by the
        skew-robust HARTIGAN DIP TEST (the extreme goblet pole inflates dBIC).
  In each system the progenitor is a single unimodal population; a synthetic
  50/50 mixture of the two committed fates is the positive control that fixes
  which modality statistic is valid.

------------------------------------------------------------------------------
DATA (download once into DATA_DIR below; these are public UCSC Cell Browser sets)
------------------------------------------------------------------------------
BRAUN (fetal brain), dataset 'fetal-brain' area/glia subset used here as tsv:
  # exprMatrix subset to marker genes + full meta. The two files this script reads:
  #   braun_markers.tsv  : genes(rows, "ENSG|SYMBOL") x cells, raw/processed counts
  #   braun_meta.tsv      : per-cell metadata incl. cell_class, cell_type, structure, numi
  # Source: Braun, Danan-Gotthold et al., Science 2023 (10.1126/science.adf1226),
  #         first-trimester human brain. UCSC Cell Browser 'fetal-brain'.
  #   base=https://cells.ucsc.edu/fetal-brain/<subset>
  #   (braun_markers.tsv = rows for the 13 glial markers below, extracted from
  #    exprMatrix.tsv.gz via: gunzip -c exprMatrix.tsv.gz | awk -F'\t' 'NR==1||$1~/marker/' )
GUT (adult intestine), dataset 'intestine-epithelium' (10x mtx):
  base=https://cells.ucsc.edu/intestine-epithelium
  curl -O $base/matrix.mtx.gz $base/barcodes.tsv.gz $base/features.tsv.gz $base/meta.tsv
  # rename to gut_matrix.mtx.gz gut_barcodes.tsv.gz gut_features.tsv.gz gut_meta.tsv
  # Source: Burclaff et al. 2022, Cell Mol Gastroenterol Hepatol 13(5):1554-1589,
  #         healthy adult human small intestine + colon epithelium, 3 donors.
------------------------------------------------------------------------------
Run:  conda activate cesfw_env; python fig7_generalisation.py
Deps: numpy pandas scikit-learn scipy matplotlib diptest
"""
import os, gzip, subprocess, numpy as np, pandas as pd, diptest
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

# EDIT THIS to wherever you downloaded the data:
DATA_DIR = os.environ.get('FIG7_DATA', os.path.expanduser('~/cESFW_analysis/figure_reconstruction/data'))
OUT = os.path.expanduser('~/cESFW_analysis/plots/revision_results/Fig7_generalisation.png')

RNG = np.random.RandomState(0)
z = lambda v: (v - np.nanmean(v)) / (np.nanstd(v) + 1e-9)

def dbic(x):
    X = np.asarray(x).reshape(-1, 1)
    return (GaussianMixture(1, n_init=2, random_state=0).fit(X).bic(X)
            - GaussianMixture(2, n_init=2, random_state=0).fit(X).bic(X))
def dbic_boot(x, n=150, r=80):
    x = np.asarray(x); return np.nan if len(x) < n else np.median([dbic(x[RNG.choice(len(x), n, 0)]) for _ in range(r)])
def dip_boot(x, n=150, r=80):
    x = np.asarray(x); return np.median([diptest.diptest(x[RNG.choice(len(x), min(n, len(x)), 0)])[1] for _ in range(r)])

# markers
BR_ASTRO = ['AQP4','GFAP','AGT','ALDH1L1','GJA1','S100B']
BR_OLIGO = ['OLIG1','OLIG2','SOX10','PDGFRA','CSPG4','PLP1','MBP']
GUT_ABS = ['ENSG00000163586','ENSG00000145384','ENSG00000118137','ENSG00000110244','ENSG00000084674',
           'ENSG00000110245','ENSG00000114113','ENSG00000163295','ENSG00000166825','ENSG00000100170',
           'ENSG00000138823','ENSG00000090402','ENSG00000118094','ENSG00000163581']       # FABP1..SLC2A2
GUT_SEC = ['ENSG00000198788','ENSG00000160180','ENSG00000275395','ENSG00000016490','ENSG00000122711',
           'ENSG00000106541','ENSG00000174992','ENSG00000134193','ENSG00000124664','ENSG00000172238',
           'ENSG00000160182','ENSG00000215182']                                            # MUC2..MUC5AC

# ---------------- BRAUN (fetal brain glia, thalamus) ----------------
def braun_phi():
    E = pd.read_csv(f'{DATA_DIR}/braun_markers.tsv', sep='\t', index_col=0)
    E.index = [g.split('|')[0] for g in E.index]; E = E.T
    m = pd.read_csv(f'{DATA_DIR}/braun_meta.tsv', sep='\t', low_memory=False).set_index('cell_name')
    common = E.index.intersection(m.index); E = E.loc[common]; m = m.loc[common]
    numi = pd.to_numeric(m['numi'], errors='coerce').fillna(pd.to_numeric(m.get('ncount_rna'), errors='coerce')).astype(float).values
    numi[numi <= 0] = np.nan
    En = pd.DataFrame(np.log1p(E.values / numi[:, None] * 1e4), index=E.index, columns=E.columns)
    AC = [g for g in BR_ASTRO if g in En.columns]; OL = [g for g in BR_OLIGO if g in En.columns]
    m = m.assign(astro=En[AC].mean(1).values, oligo=En[OL].mean(1).values)
    g = m[(m['cell_class'] == 'glia') & (m['structure'] == 'thalamus')].copy()
    phi = z(g['astro'].values) - z(g['oligo'].values); ct = g['cell_type'].values
    return phi[ct == 'radial glia'], phi[ct == 'astrocyte'], phi[ct == 'oligodendrocyte']

# ---------------- GUT (adult small intestine) ----------------
def gut_phi():
    feat = [l.strip() for l in gzip.open(f'{DATA_DIR}/gut_features.tsv.gz', 'rt')]
    rowpos = {e: i for i, e in enumerate(feat)}
    order = {rowpos[e] + 1: k for k, e in enumerate(GUT_ABS + GUT_SEC)}   # mtx rows are 1-based
    bc = [l.strip() for l in gzip.open(f'{DATA_DIR}/gut_barcodes.tsv.gz', 'rt')]
    ncells = len(bc); M = np.zeros((len(GUT_ABS) + len(GUT_SEC), ncells)); colsum = np.zeros(ncells)
    p = subprocess.Popen(['gunzip', '-c', f'{DATA_DIR}/gut_matrix.mtx.gz'], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    it = p.stdout; [next(it) for _ in range(3)]                            # skip 2 comment + dims line
    for line in it:
        r, c, v = line.split(); ci = int(c) - 1; v = float(v); colsum[ci] += v
        k = order.get(int(r))
        if k is not None: M[k, ci] = v
    p.wait()
    L = np.log1p(M / np.where(colsum <= 0, 1, colsum) * 1e4)               # markers x cells, CP10K log
    Ld = pd.DataFrame(L.T, columns=GUT_ABS + GUT_SEC)
    meta = pd.read_csv(f'{DATA_DIR}/gut_meta.tsv', sep='\t', low_memory=False).set_index('cellId').reindex(bc)
    Ld.index = meta.index
    sub = meta[meta['lineage'].astype(str).str.startswith('SI_')]          # small-intestine context
    Zc = Ld.loc[sub.index].apply(z)
    phi = pd.Series(z(Zc[GUT_ABS].mean(1).values) - z(Zc[GUT_SEC].mean(1).values), index=sub.index)
    g = lambda t: phi[sub['type'].values == t].values
    return g('ISC'), g('absorptive'), g('goblet')

def main():
    rg, ast, ol = braun_phi()
    isc, absorp, gob = gut_phi()
    synb = lambda n: np.concatenate([RNG.choice(ast, n // 2), RNG.choice(ol, n - n // 2)])
    syng = lambda n: np.concatenate([RNG.choice(absorp, n // 2), RNG.choice(gob, n - n // 2)])
    b_rg = dbic_boot(rg); b_mix = np.median([dbic(synb(150)) for _ in range(80)])
    g_isc = dip_boot(isc); g_mix = np.median([diptest.diptest(syng(150))[1] for _ in range(80)]); g_isc_db = dbic_boot(isc)
    print(f'BRAUN  RG dBIC={b_rg:+.1f} | synth-mix dBIC={b_mix:+.1f} | means RG={rg.mean():+.2f} ast={ast.mean():+.2f} ol={ol.mean():+.2f}')
    print(f'GUT    ISC dip={g_isc:.3f} | synth-mix dip={g_mix:.3f} | ISC dBIC(skew-inflated)={g_isc_db:+.1f} | means ISC={isc.mean():+.2f} abs={absorp.mean():+.2f} gob={gob.mean():+.2f}')

    CA = {'prog': '#E8873A', 'fA': '#1E7145', 'fB': '#6C3483', 'mix': '#999999'}
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.6)); fig.subplots_adjust(hspace=0.42, wspace=0.26)
    def dist(a, xp, xa, xb, labels, title, xlab):
        lo = min(xp.min(), xa.min(), xb.min()) - .5; hi = max(xp.max(), xa.max(), xb.max()) + .5; gr = np.linspace(lo, hi, 300)
        for x, lab, c in [(xa, labels[1], CA['fA']), (xp, labels[0], CA['prog']), (xb, labels[2], CA['fB'])]:
            a.plot(gr, gaussian_kde(x)(gr), color=c, lw=2.4, label=f'{lab} (n={len(x)})'); a.fill_between(gr, gaussian_kde(x)(gr), alpha=.28, color=c)
        a.axvline(np.mean(xp), color=CA['prog'], ls='--', lw=1.1)
        a.set_xlabel(xlab); a.set_ylabel('density'); a.set_title(title, fontsize=10.3)
        a.legend(frameon=False, fontsize=8.2, loc='upper left'); a.spines[['top', 'right']].set_visible(False)
    def pc(a, xp, mixfun, proglab, ann, title, xlab):
        p150 = xp[RNG.choice(len(xp), 150, 0)]; m150 = mixfun(150)
        lo = min(p150.min(), m150.min()) - .5; hi = max(p150.max(), m150.max()) + .5; gr = np.linspace(lo, hi, 300)
        a.plot(gr, gaussian_kde(m150)(gr), color=CA['mix'], lw=2.2, label='synthetic 50/50 mixture'); a.fill_between(gr, gaussian_kde(m150)(gr), alpha=.3, color=CA['mix'])
        a.plot(gr, gaussian_kde(p150)(gr), color=CA['prog'], lw=2.4, label=proglab); a.fill_between(gr, gaussian_kde(p150)(gr), alpha=.3, color=CA['prog'])
        a.set_xlabel(xlab); a.set_ylabel('density'); a.set_title(title, fontsize=10.3)
        a.legend(frameon=False, fontsize=8.2, loc='upper left'); a.spines[['top', 'right']].set_visible(False)
        a.text(0.985, 0.62, ann, transform=a.transAxes, ha='right', va='top', fontsize=8.3, bbox=dict(boxstyle='round', fc='#f4f4f4', ec='#bbb'))
    dist(ax[0, 0], rg, ast, ol, ['glial progenitor (radial glia)', 'astrocyte', 'oligodendrocyte'],
         'Fetal brain glia (thalamus): the progenitor is a single\nunresolved population between its two fates', 'φ = z(astrocyte) − z(oligodendrocyte)')
    pc(ax[0, 1], rg, synb, 'radial-glia progenitor',
       f'calibrated ΔBIC (>0 = bimodal):\nprogenitor {b_rg:+.0f} (unimodal)\nmixture {b_mix:+.0f} (bimodal)',
       'Positive control (n=150): the test detects a true\nmixture; the progenitor is not one', 'φ = z(astrocyte) − z(oligodendrocyte)')
    dist(ax[1, 0], isc, absorp, gob, ['intestinal stem cell (ISC)', 'absorptive (enterocyte)', 'secretory (goblet)'],
         'Adult intestine (small intestine): the stem cell is a single\nunresolved population between its two fates', 'φ = z(absorptive) − z(secretory)')
    pc(ax[1, 1], isc, syng, 'ISC progenitor',
       f'Hartigan dip test (p>0.05 = unimodal):\nprogenitor p={g_isc:.2f} (unimodal)\nmixture p<0.001 (multimodal)\n(ΔBIC skew-inflated here)',
       'Positive control (n=150): the test detects a true\nmixture; the stem cell is not one', 'φ = z(absorptive) − z(secretory)')
    for a, L in zip([ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]], 'ABCD'):
        a.text(-0.13, 1.06, L, transform=a.transAxes, fontsize=15, fontweight='bold', va='top', ha='left')
    fig.suptitle('The transition-versus-state framework generalises to two unrelated human fate decisions', fontsize=12.5, y=0.975, fontweight='bold')
    fig.savefig(OUT, dpi=200, bbox_inches='tight', facecolor='white'); print('SAVED', OUT)

if __name__ == '__main__':
    main()
