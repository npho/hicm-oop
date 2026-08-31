#!/usr/bin/env python
# =============================================================================
#  Figure S14  --  *** RECONSTRUCTION, NOT ORIGINAL CODE ***
# =============================================================================
#  The script that produced Figure S14 (`Fig_negative_control.png`, written
#  2026-07-19 09:25) does NOT exist on disk. Exhaustive search: no .py, .ipynb,
#  .sh, .txt or .md file anywhere under cESFW_analysis/ or cesfw_runs/
#  references that filename, and no Python file was modified on that date other
#  than the unrelated generalisation tests (lung/trophoblast/cardiac/gbm).
#  It was evidently run ad hoc and not saved.
#
#  An earlier version of this audit package attributed the figure to
#  `figure_negative_controls.py`. THAT WAS WRONG - that script builds naive-PSC
#  projection controls (t2iLGo/PXGL/Mazid/Bi/EPSC) and writes
#  `negative_controls_extended.png`. It has been removed from this folder.
#
#  This file re-implements the analysis exactly as the published caption and
#  the Methods describe it, so the auditor can check the published numbers.
#  It is a reconstruction: it reproduces the RESULT, and should not be cited as
#  the code that made the figure.
#
#  PUBLISHED CAPTION (targets to reproduce):
#    "(A) Hartigan dip test p for real populations: single cell types, the ICM,
#     epiblast, hypoblast and trophectoderm, are each unimodal (p>0.05), whereas
#     a population formed by pooling real committed epiblast and hypoblast cells
#     is bimodal (p<0.05) ... (B) the ICM (dip p=0.95) is a single population,
#     not a resolved epiblast+hypoblast mixture (dip p=0.02). DeltaBIC is
#     uninformative at large n (trophectoderm, n=993, DeltaBIC>1000 yet dip
#     p=0.998)"
#
#  Method (manuscript Methods, cESFW-independence section): standard total-count
#  normalisation to 10,000 + log1p over the FULL transcriptome (deliberately NOT
#  the cESFW gene set), independent Stirparo et al. 2018 lineage labels, and
#  phi = z(epiblast score) - z(hypoblast score) from canonical markers.
#
#  Usage:  python RECONSTRUCTED_negative_control.py
# =============================================================================
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
sys.path.insert(0, '/Users/adelosan/cESFW_analysis')  # repoint if relocating
import _load_data_yanagida as ldy

try:
    import diptest
except ImportError:
    sys.exit('diptest required:  pip install diptest')
from sklearn.mixture import GaussianMixture

EPI = ['NANOG', 'SOX2', 'KLF17', 'TDGF1', 'DPPA5', 'GDF3', 'PRDM14']
HYP = ['GATA6', 'GATA4', 'SOX17', 'PDGFRA', 'FOXA2', 'HNF1B', 'APOA1', 'COL4A1']
SEED = 0
N_SUB = 150          # matches the subsampling used elsewhere in the paper
N_BOOT = 80

Sample_Info, Counts, _, _, _, _ = ldy.load_embryo_data()

# standard normalisation over the full transcriptome (NOT cESFW-restricted)
lib = Counts.sum(axis=1)
lib[lib == 0] = 1
lc = np.log1p(Counts.div(lib, axis=0) * 1e4)

epi = [g for g in EPI if g in lc.columns]
hyp = [g for g in HYP if g in lc.columns]
print(f'epiblast markers used  ({len(epi)}/{len(EPI)}): {epi}')
print(f'hypoblast markers used ({len(hyp)}/{len(HYP)}): {hyp}\n')


def z(x):
    s = x.std()
    return (x - x.mean()) / (s if s else 1.0)


phi = z(lc[epi].mean(axis=1)) - z(lc[hyp].mean(axis=1))

lab = Sample_Info['Stirparo_Labels']
print('Stirparo label counts:')
print(lab.value_counts().to_string(), '\n')

GROUPS = {
    'ICM': ['ICM', 'Early ICM'],
    'Epiblast': ['Epiblast'],
    'Hypoblast': ['Hypoblast'],
    'Trophectoderm': ['TE'],
}


def dip_p(v, sub=N_SUB, boot=N_BOOT, seed=SEED):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) <= sub:
        return diptest.diptest(v)[1], len(v)
    rng = np.random.default_rng(seed)
    ps = [diptest.diptest(rng.choice(v, sub, replace=False))[1] for _ in range(boot)]
    return float(np.median(ps)), len(v)


def dbic(v):
    v = np.asarray(v, float).reshape(-1, 1)
    v = v[np.isfinite(v[:, 0])]
    if len(v) < 20:
        return np.nan
    b1 = GaussianMixture(1, n_init=3, random_state=SEED).fit(v).bic(v)
    b2 = GaussianMixture(2, n_init=3, random_state=SEED).fit(v).bic(v)
    return b1 - b2


print('=' * 70)
print('PANEL A - real populations, each tested for unimodality')
print('=' * 70)
rows = []
for name, labs in GROUPS.items():
    m = lab.isin(labs).values
    p, n = dip_p(phi[m])
    rows.append({'population': name, 'n': n, 'dip_p': round(p, 3),
                 'dBIC': round(dbic(phi[m]), 1),
                 'verdict': 'unimodal' if p > 0.05 else 'BIMODAL'})

# the negative control: pool REAL committed epiblast + hypoblast cells
mix = lab.isin(['Epiblast', 'Hypoblast']).values
p, n = dip_p(phi[mix])
rows.append({'population': 'POOLED Epi+Hyp (control)', 'n': n, 'dip_p': round(p, 3),
             'dBIC': round(dbic(phi[mix]), 1),
             'verdict': 'unimodal' if p > 0.05 else 'BIMODAL'})
print(pd.DataFrame(rows).to_string(index=False))

print('\n' + '=' * 70)
print('CHECK AGAINST THE PUBLISHED CAPTION')
print('=' * 70)
got = {r['population']: r for r in rows}
checks = [
    ('ICM unimodal, dip p ~ 0.95', got['ICM']['dip_p'], 0.95, 0.15),
    ('pooled Epi+Hyp bimodal, dip p ~ 0.02',
     got['POOLED Epi+Hyp (control)']['dip_p'], 0.02, 0.05),
    ('trophectoderm n = 993', got['Trophectoderm']['n'], 993, 0),
    ('trophectoderm dip p ~ 0.998', got['Trophectoderm']['dip_p'], 0.998, 0.05),
]
for label, obs, exp, tol in checks:
    ok = abs(obs - exp) <= tol
    print(f'  [{"MATCH " if ok else "DIFFER"}] {label:42s} observed={obs}')
tro = got['Trophectoderm']
print(f'  [{"MATCH " if tro["dBIC"] > 1000 else "DIFFER"}] trophectoderm dBIC > 1000 '
      f'(dip still says unimodal)      observed={tro["dBIC"]}')
print('\nA full set of MATCH lines means this reconstruction reproduces the '
      'published figure\'s numbers.\nIt remains a reconstruction, not the '
      'original code.')
