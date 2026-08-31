# ============================================================================
# Figure S15  --  homeostatic lung (Sikkema 2023)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The KRT8+ AT2-to-AT1 alveolar intermediate is a single transition in
#   both homeostatic and fibrotic lung.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   homeostatic lung (Sikkema 2023)
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S15. The KRT8-positive AT2-to-AT1 alveolar intermediate is a
#   single transition in both homeostatic and fibrotic lung. Order-parameter
#   modality test (φ = z(AT1) − z(AT2)) applied to (A) the healthy Human
#   Lung Cell Atlas (Sikkema et al., 2023) and (B) a pulmonary-fibrosis
#   atlas (Habermann et al., 2020). In each, the transitional intermediate
#   (homeostatic 'AT0'; fibrotic 'Transitional AT2') is a single unimodal
#   population between the AT2 and AT1 poles (dip p=0.94 and 0.80), while a
#   synthetic AT2+AT1 mixture is bimodal (dip p<0.001), confirming the test
#   discriminates. The debated KRT8-positive intermediate is therefore a
#   transient AT2-to-AT1 transition, not a discrete cell state, in both
#   homeostasis and the fibrotic context where the question is most
#   contested.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/lung_at2at1_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""
Adjudicate: is the KRT8+ alveolar regeneration intermediate (transitional AT2 /
ADI / PATS / AT0) a DISCRETE cell state or a single unresolved TRANSITION along
the AT2 -> AT1 axis?

REGENERATION / 1->1 conversion adaptation of the bipotent-transition modality test.
Poles = AT2 and AT1.  Transitional = HLCA "AT0" annotation.
Question: do transitional cells form ONE unimodal distribution INTERMEDIATE
between the two poles (-> TRANSITION), or their own separated mode / a resolved
AT2+AT1 mixture (-> DISCRETE STATE)?

Dataset: The integrated Human Lung Cell Atlas (HLCA), core (healthy reference),
Sikkema et al. 2023 Nat Med, 10.1038/s41591-023-02327-2. CELLxGENE Discover
dataset 066943a2-fdac-4b29-b348-40cede398e4e.

phi = z(AT1 markers) - z(AT2 markers), z-scored across AT2+AT0+AT1 cells.
"""
import numpy as np, pandas as pd, anndata as ad
import diptest
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.RandomState(0)
H5 = "/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad/lung_hlca_core.h5ad"
OUTFIG = "/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_lung_at2at1.png"

AT2M = ['SFTPC','SFTPB','SFTPA1','SFTPA2','LAMP3','ETV5','NAPSA']
AT1M = ['AGER','PDPN','HOPX','CAV1','RTKN2','CLIC5','EMP2']
VERIFY = ['KRT8','CLDN4','SFTPC','SFTPB','AGER','PDPN','HOPX']

z = lambda v: (v - np.mean(v)) / (np.std(v) + 1e-9)

def dbic1(x):
    X = np.asarray(x, float).reshape(-1,1)
    return (GaussianMixture(1, n_init=2, random_state=0).fit(X).bic(X)
            - GaussianMixture(2, n_init=2, random_state=0).fit(X).bic(X))

def boot(x, n=150, reps=80):
    """calibrated (fixed n) median dBIC, %bimodal, median dip-p."""
    x = np.asarray(x, float)
    if len(x) < n:
        return np.nan, np.nan, np.nan
    db, dp = [], []
    for _ in range(reps):
        s = x[rng.choice(len(x), n, replace=False)]
        db.append(dbic1(s))
        dp.append(diptest.diptest(s)[1])
    db = np.array(db)
    return np.median(db), np.mean(db > 0), np.median(dp)

def gmm2(x):
    X = np.asarray(x, float).reshape(-1,1)
    g = GaussianMixture(2, n_init=5, random_state=0).fit(X)
    o = np.argsort(g.means_.ravel())
    mu = g.means_.ravel()[o]; w = g.weights_[o]
    sd = np.sqrt(g.covariances_.ravel()[o])
    sep = abs(mu[1]-mu[0]) / np.sqrt((sd[0]**2 + sd[1]**2)/2)
    return mu, w, sd, sep

# ---------------------------------------------------------------- load
print("loading HLCA core (backed)...")
A = ad.read_h5ad(H5, backed='r')
ann = A.obs['ann_finest_level'].astype(str).values
groups = {'AT2':'AT2', 'AT0':'AT0', 'AT1':'AT1'}   # pole/trans/pole
mask = np.isin(ann, list(groups.keys()))
print("selected cells:", mask.sum(), "of", A.n_obs)

# gene index via var feature_name
fn = A.var['feature_name'].astype(str).values
gidx = {g:i for i,g in enumerate(fn)}
need = AT2M + AT1M + [g for g in VERIFY if g not in AT2M+AT1M]
need = [g for g in need if g in gidx]
cols = [gidx[g] for g in need]

# raw counts for selected cells (raw slot = counts per CELLxGENE schema)
raw = A.raw
sub_counts = raw.X[np.where(mask)[0], :]        # cells x all genes (sparse)
tot = np.asarray(sub_counts.sum(1)).ravel()     # total counts per cell (all genes)
mat = np.asarray(sub_counts[:, cols].todense())  # cells x needed genes (raw)
del sub_counts
# CP10K + log1p
cp10k = mat / tot[:,None] * 1e4
logn = np.log1p(cp10k)
expr = pd.DataFrame(logn, columns=need)
lab = ann[mask]

# ---------------------------------------------------------------- phi
at2_present = [g for g in AT2M if g in expr.columns]
at1_present = [g for g in AT1M if g in expr.columns]
at2score = expr[at2_present].mean(1).values
at1score = expr[at1_present].mean(1).values
phi = z(at1score) - z(at2score)   # z across all AT2+AT0+AT1 cells

df = pd.DataFrame({'lab':lab, 'phi':phi})
for g in VERIFY:
    if g in expr.columns: df[g] = expr[g].values

# ---------------------------------------------------------------- verify transitional identity
print("\n=== transitional-identity verification (mean log-norm expression) ===")
vt = df.groupby('lab')[[g for g in VERIFY if g in df.columns]].mean().reindex(['AT2','AT0','AT1'])
print(vt.round(2).to_string())
print("\nmean phi by group:")
print(df.groupby('lab')['phi'].agg(['mean','median','count']).reindex(['AT2','AT0','AT1']).round(3).to_string())

# ---------------------------------------------------------------- positive control: synthetic 50/50 AT2+AT1
nmix = min((lab=='AT1').sum(), (lab=='AT2').sum(), 4000)
a2 = rng.choice(np.where(lab=='AT2')[0], nmix, replace=False)
a1 = rng.choice(np.where(lab=='AT1')[0], nmix, replace=False)
mix = phi[np.concatenate([a2,a1])]

# ---------------------------------------------------------------- modality table
print("\n=== calibrated modality (subsample n=150, 80 reps) ===")
rows = []; B = {}
for name, x in [('AT2 pole', phi[lab=='AT2']),
                ('AT0 transitional', phi[lab=='AT0']),
                ('AT1 pole', phi[lab=='AT1']),
                ('synthetic 50/50 mix', mix)]:
    md, fb, dp = boot(x)
    B[name] = (md, fb, dp)
    verdict = 'BIMODAL' if (md>0) else 'unimodal'
    rows.append((name, len(x), np.mean(x), np.median(x), dp, md, fb, verdict))
    print(f"{name:22s} n={len(x):6d} meanphi={np.mean(x):+.2f} dip_p={dp:.3f} "
          f"dBIC={md:+.1f} ({100*fb:.0f}% bimodal) -> {verdict}")
tab = pd.DataFrame(rows, columns=['pop','n','mean_phi','median_phi','dip_p','dBIC_median','frac_bimodal','verdict'])
tab.to_csv("/Users/adelosan/cESFW_analysis/plots/revision_results/lung_at2at1_table.csv", index=False)

# ---------------------------------------------------------------- transitional 2-comp GMM decomposition
mu, w, sd, sep = gmm2(phi[lab=='AT0'])
print("\n=== AT0 transitional 2-component GMM ===")
print(f"comp1 mean={mu[0]:+.2f} sd={sd[0]:.2f} w={w[0]:.2f}")
print(f"comp2 mean={mu[1]:+.2f} sd={sd[1]:.2f} w={w[1]:.2f}")
print(f"standardized separation = {sep:.2f}  (pole means: AT2={df[df.lab=='AT2'].phi.mean():+.2f}, AT1={df[df.lab=='AT1'].phi.mean():+.2f})")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(13,5.2))
cols_g = {'AT2':'#2c7fb8','AT0':'#d95f0e','AT1':'#31a354'}
bins = np.linspace(np.percentile(phi,0.3), np.percentile(phi,99.7), 60)
for g in ['AT2','AT0','AT1']:
    ax[0].hist(phi[lab==g], bins=bins, density=True, alpha=0.55,
               color=cols_g[g], label=f'{g} (n={(lab==g).sum()})')
for g in ['AT2','AT0','AT1']:
    ax[0].axvline(df[df.lab==g].phi.mean(), color=cols_g[g], ls='--', lw=1.5)
ax[0].set_xlabel(r'$\phi$  =  z(AT1 program) - z(AT2 program)')
ax[0].set_ylabel('density')
ax[0].legend(frameon=False, fontsize=9)
ax[0].set_title('A. Transitional (AT0) sits between the AT2 and AT1 poles', fontsize=11, loc='left')

md_t, fb_t, dp_t = B['AT0 transitional']
md_m, fb_m, dp_m = B['synthetic 50/50 mix']
ax[1].hist(phi[lab=='AT0'], bins=bins, density=True, alpha=0.6, color=cols_g['AT0'],
           label=f'AT0 transitional\ndip p={dp_t:.2f}, dBIC={md_t:+.1f}')
ax[1].hist(mix, bins=bins, density=True, alpha=0.45, color='#555555',
           label=f'synthetic 50/50 AT2+AT1\ndip p={dp_m:.3f}, dBIC={md_m:+.1f}')
ax[1].set_xlabel(r'$\phi$')
ax[1].set_ylabel('density')
ax[1].legend(frameon=False, fontsize=9)
verdict_t = 'unimodal (a transition)' if md_t<0 else 'bimodal (a mixture)'
ax[1].set_title(f'B. Transitional is {verdict_t}; positive control is bimodal', fontsize=11, loc='left')
plt.tight_layout()
plt.savefig(OUTFIG, dpi=150)
print("\nsaved figure:", OUTFIG)
