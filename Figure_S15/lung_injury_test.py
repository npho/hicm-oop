# ============================================================================
# Figure S15  --  fibrotic lung (Habermann 2020)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The KRT8+ AT2-to-AT1 alveolar intermediate is a single transition in
#   both homeostatic and fibrotic lung.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   fibrotic lung (Habermann 2020)
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/lung_injury_test.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""
INJURY / FIBROSIS confirmation of the KRT8+ AT2->AT1 transitional-intermediate test.

Companion to the HEALTHY Human Lung Cell Atlas analysis (lung_at2at1_test.py),
which showed the KRT8+/CLDN4+ "AT0" AT2->AT1 intermediate is a single unresolved
TRANSITION (unimodal phi, intermediate between AT2 and AT1 poles; synthetic
AT2+AT1 mixture bimodal positive control). Here we ask the same question where
the debate is loudest -- in DISEASED / FIBROTIC lung.

Question: does the DISEASE-associated KRT8+ AT2->AT1 intermediate read as a
single TRANSITION (unimodal, intermediate) or a DISCRETE STATE (own separated
mode / resolves to the poles)?

Dataset : Habermann et al. 2020, Sci Adv 6:eaba1972 (GSE135893), "Single-cell
          RNA-sequencing reveals profibrotic roles of distinct epithelial and
          mesenchymal lineages in pulmonary fibrosis". Human distal lung,
          IPF/fibrotic-ILD vs control. Author annotations include an explicit
          "Transitional AT2" (KRT8-high) cluster, distinct from the fibrosis-
          specific aberrant-basaloid "KRT5-/KRT17+" cluster (EXCLUDED).
Cells   : DISEASE (Status=='Disease') AT2 / Transitional AT2 / AT1 only.
Poles   : AT2 and AT1.  Transitional = author "Transitional AT2".
phi     : z(AT1 markers) - z(AT2 markers), z-scored across AT2+Transitional+AT1.

Pipeline identical to companion: CP10K + log1p over all genes; calibrated
modality (subsample n=150, 80 boots) reporting median dBIC, %bimodal, median
Hartigan dip p; 2-component GMM decomposition of the transitional population;
synthetic 50/50 AT2+AT1 mixture positive control.
"""
import numpy as np, pandas as pd
import diptest
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.RandomState(0)
D   = "/private/tmp/claude-501/-Users-adelosan/73724e25-f2ba-48fa-a729-bbdde91c906b/scratchpad/"
NPZ = D + "lunginj_extracted.npz"
OUTFIG = "/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_lung_injury.png"
OUTTAB = "/Users/adelosan/cESFW_analysis/plots/revision_results/lung_injury_table.csv"

AT2M = ['SFTPC','SFTPB','SFTPA1','SFTPA2','LAMP3','ETV5','NAPSA']
AT1M = ['AGER','PDPN','HOPX','CAV1','RTKN2','CLIC5','EMP2']
VERIFY = ['KRT8','CLDN4','SFTPC','SFTPB','AGER','PDPN','HOPX']

z = lambda v: (v - np.mean(v)) / (np.std(v) + 1e-9)

def dbic1(x):
    X = np.asarray(x, float).reshape(-1,1)
    return (GaussianMixture(1, n_init=2, random_state=0).fit(X).bic(X)
            - GaussianMixture(2, n_init=2, random_state=0).fit(X).bic(X))

def boot(x, n=150, reps=80):
    x = np.asarray(x, float)
    if len(x) < n:
        # fall back to full-sample single estimate if pole is small
        return dbic1(x), float(dbic1(x) > 0), diptest.diptest(x)[1]
    db, dp = [], []
    for _ in range(reps):
        s = x[rng.choice(len(x), n, replace=False)]
        db.append(dbic1(s)); dp.append(diptest.diptest(s)[1])
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

# --------------------------------------------------------------- load extracted
Z = np.load(NPZ, allow_pickle=True)
mat = Z['mat']; tot = Z['tot']; lab = Z['lab'].astype(str); genes = list(Z['genes'].astype(str))
# CP10K + log1p (total = counts over ALL genes, computed during extraction)
logn = np.log1p(mat / tot[:,None] * 1e4)
expr = pd.DataFrame(logn, columns=genes)

at2_present = [g for g in AT2M if g in expr.columns]
at1_present = [g for g in AT1M if g in expr.columns]
at2score = expr[at2_present].mean(1).values
at1score = expr[at1_present].mean(1).values
phi = z(at1score) - z(at2score)

df = pd.DataFrame({'lab':lab, 'phi':phi})
for g in VERIFY:
    if g in expr.columns: df[g] = expr[g].values

LABELS = ['AT2','Transitional AT2','AT1']

print("AT2 markers used:", at2_present)
print("AT1 markers used:", at1_present)
print("\n=== transitional-identity verification (mean log-norm expression) ===")
vt = df.groupby('lab')[[g for g in VERIFY if g in df.columns]].mean().reindex(LABELS)
print(vt.round(2).to_string())
print("\nmean phi by group:")
print(df.groupby('lab')['phi'].agg(['mean','median','count']).reindex(LABELS).round(3).to_string())

# --------------------------------------------------------------- positive control
nmix = min((lab=='AT1').sum(), (lab=='AT2').sum())
a2 = rng.choice(np.where(lab=='AT2')[0], nmix, replace=False)
a1 = rng.choice(np.where(lab=='AT1')[0], nmix, replace=False)
mix = phi[np.concatenate([a2,a1])]

print("\n=== calibrated modality (subsample n=150, 80 reps) ===")
rows = []; B = {}
for name, x in [('AT2 pole', phi[lab=='AT2']),
                ('Transitional AT2', phi[lab=='Transitional AT2']),
                ('AT1 pole', phi[lab=='AT1']),
                ('synthetic 50/50 mix', mix)]:
    md, fb, dp = boot(x)
    B[name] = (md, fb, dp)
    verdict = 'BIMODAL' if (md>0) else 'unimodal'
    rows.append((name, len(x), np.mean(x), np.median(x), dp, md, fb, verdict))
    print(f"{name:22s} n={len(x):6d} meanphi={np.mean(x):+.2f} dip_p={dp:.3f} "
          f"dBIC={md:+.1f} ({100*fb:.0f}% bimodal) -> {verdict}")
pd.DataFrame(rows, columns=['pop','n','mean_phi','median_phi','dip_p','dBIC_median',
                            'frac_bimodal','verdict']).to_csv(OUTTAB, index=False)

# --------------------------------------------------------------- transitional GMM
mu, w, sd, sep = gmm2(phi[lab=='Transitional AT2'])
a2m = df[df.lab=='AT2'].phi.mean(); a1m = df[df.lab=='AT1'].phi.mean()
print("\n=== Transitional AT2  2-component GMM ===")
print(f"comp1 mean={mu[0]:+.2f} sd={sd[0]:.2f} w={w[0]:.2f}")
print(f"comp2 mean={mu[1]:+.2f} sd={sd[1]:.2f} w={w[1]:.2f}")
print(f"standardized separation = {sep:.2f}  (pole means: AT2={a2m:+.2f}, AT1={a1m:+.2f})")
inzone = (min(mu)>a2m) and (max(mu)<a1m)
print("both GMM comps strictly between the poles:", inzone)

# --------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 2, figsize=(13,5.2))
cg = {'AT2':'#2c7fb8','Transitional AT2':'#d95f0e','AT1':'#31a354'}
bins = np.linspace(np.percentile(phi,0.3), np.percentile(phi,99.7), 55)
for g in LABELS:
    ax[0].hist(phi[lab==g], bins=bins, density=True, alpha=0.55, color=cg[g],
               label=f'{g} (n={(lab==g).sum()})')
for g in LABELS:
    ax[0].axvline(df[df.lab==g].phi.mean(), color=cg[g], ls='--', lw=1.5)
ax[0].set_xlabel(r'$\phi$  =  z(AT1 program) - z(AT2 program)')
ax[0].set_ylabel('density'); ax[0].legend(frameon=False, fontsize=9)
ax[0].set_title('A. IPF/fibrotic lung: Transitional AT2 sits between AT2 and AT1 poles',
                fontsize=10.5, loc='left')

md_t, fb_t, dp_t = B['Transitional AT2']
md_m, fb_m, dp_m = B['synthetic 50/50 mix']
ax[1].hist(phi[lab=='Transitional AT2'], bins=bins, density=True, alpha=0.6, color=cg['Transitional AT2'],
           label=f'Transitional AT2 (disease)\ndip p={dp_t:.2f}, dBIC={md_t:+.1f}')
ax[1].hist(mix, bins=bins, density=True, alpha=0.45, color='#555555',
           label=f'synthetic 50/50 AT2+AT1\ndip p={dp_m:.3f}, dBIC={md_m:+.1f}')
ax[1].set_xlabel(r'$\phi$'); ax[1].set_ylabel('density'); ax[1].legend(frameon=False, fontsize=9)
verdict_t = 'unimodal (a transition)' if md_t<0 else 'bimodal (a mixture)'
ax[1].set_title(f'B. Transitional AT2 is {verdict_t}; positive control is bimodal',
                fontsize=10.5, loc='left')
plt.tight_layout(); plt.savefig(OUTFIG, dpi=150)
print("\nsaved figure:", OUTFIG)
print("saved table :", OUTTAB)
