# ============================================================================
# Figure 6  --  marmoset ICM bifurcation
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mouse ICM is already resolved into two populations; the human and
#   marmoset ICM are not. The species contrast that anchors the transition
#   claim.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   marmoset ICM bifurcation
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 6. ICM resolution across species and along the human lineage. (A)
#   Mouse and primate ICMs differ in lineage resolution. Distributions of
#   the ICM lineage order parameter under two metrics are shown: the
#   dropout-insensitive multigene epiblast-minus-hypoblast composite (φ₃;
#   top) and the canonical NANOG–GATA6 toggle (φ₁; bottom). Kernel densities
#   are overlaid on histograms; ΔBIC summarizes support for a two-component
#   Gaussian mixture (positive = bimodal/resolved). Under φ₃, human (cESFW
#   atlas) and marmoset (Boroviak et al., 2018; ICM_CS3) ICMs are decisively
#   unimodal, whereas mouse (Nowotschin et al., 2019; E3.5 inner cells) is
#   strongly bimodal, consistent with salt-and-pepper lineage resolution.
#   The human φ₁ result lies near the unimodal–bimodal boundary (ΔBIC ≈ −3),
#   reflecting NANOG dropout (~16% of cells), whereas marmoset and mouse
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/marmoset_bifurcation.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""MARMOSET critical-signature test (primate-conservation leg of the triad).
Boroviak 2018 (E-MTAB-7078) RawCounts, de novo marker annotation (no labels in repo).
Isolate ICM/epiblast/hypoblast compartment (exclude TE + post-implantation),
phi = z(NANOG)-z(GATA6). Q: is the marmoset ICM near-critical (unimodal at
variance peak, like human) rather than switch-like (like mouse)?"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
np.seterr(all="ignore")
D = "data/nhp/marmoset/RawCounts.csv"

# row1=cellIDs(header), row2='Location ID', row3='Lineage', row4+=gene counts
meta = pd.read_csv(D, index_col=0, nrows=2)                            # rows: Location ID, Lineage
lineage = meta.iloc[1]                                                 # Lineage row: cell -> lineage
counts = pd.read_csv(D, skiprows=[1,2], index_col=0)                   # genes x cells (numeric)
counts = counts[~counts.index.duplicated()]
# preimplantation blastocyst bifurcation: ICM_CS3 -> {Epi_CS3, Hyp_CS3}
KEEP = ["ICM_CS3","Epi_CS3","Hyp_CS3"]
cells = [c for c in counts.columns if lineage.get(c) in KEEP]
ad = sc.AnnData(X=counts[cells].T.values.astype(float)); ad.var_names=counts.index; ad.obs_names=cells
ad.obs["lin"] = [lineage[c] for c in cells]
sc.pp.filter_genes(ad, min_cells=3)
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
als = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)
sub = ad
from collections import Counter
print(f"marmoset blastocyst cells: {sub.n_obs} | " + str(dict(Counter(list(ad.obs['lin'])))))

z = lambda v:(v-v.mean())/(v.std()+1e-9)
phi = z(als["NANOG"].values) - z(als["GATA6"].values)
sc.pp.neighbors(sub, n_neighbors=min(15,sub.n_obs-1), use_rep="X", metric="correlation"); sc.tl.diffmap(sub)
# root: an ICM_CS3 cell (most co-expressing / balanced)
isicm = sub.obs["lin"].values == "ICM_CS3"
coexp = z(als["NANOG"].values)+z(als["GATA6"].values)
score_root = coexp - np.abs(phi); score_root[~isicm] = -1e9
sub.uns["iroot"] = int(np.argmax(score_root)); sc.tl.dpt(sub)
pt = sub.obs["dpt_pseudotime"].values

def fita(x):
    if len(x)<20: return np.nan
    g=np.linspace(np.percentile(x,2)-.3,np.percentile(x,98)+.3,120); d=gaussian_kde(x)(g); U=-np.log(d+1e-6)
    w=d/d.sum(); Xb=np.vstack([np.ones_like(g),g,g**2,g**4]).T
    return 2*np.linalg.lstsq(Xb*w[:,None],U*w,rcond=None)[0][2]
def nmodes(x,grid):
    if len(x)<12: return 0
    d=gaussian_kde(x)(grid); return sum(1 for i in range(1,len(d)-1) if d[i]>d[i-1] and d[i]>d[i+1] and d[i]>0.05*d.max())
def bic(x):
    if len(x)<20: return np.nan
    X=x.reshape(-1,1); return GaussianMixture(1,n_init=2,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=2,random_state=0).fit(X).bic(X)

nb=6; edges=np.quantile(pt,np.linspace(0,1,nb+1)); edges[-1]+=1e-9; grid=np.linspace(phi.min()-.5,phi.max()+.5,200)
rows=[]
for i in range(nb):
    sel=(pt>=edges[i])&(pt<edges[i+1]); x=phi[sel]
    rows.append(dict(mid=pt[sel].mean(),nn=int(sel.sum()),var=x.var(),a=fita(x),nmodes=nmodes(x,grid),bic=bic(x)))
df=pd.DataFrame(rows); pk=int(np.nanargmax(df["var"].values))
print("\n=== MARMOSET per pseudotime bin ==="); print(df.round(3).to_string(index=False))
print(f"\nvariance peak: bin {pk} (pt={df.mid[pk]:.3f}); a at peak={df.a[pk]:.3f}; "
      f"min a={np.nanmin(df.a):.3f}" if np.isfinite(df.a.values).any() else "a unavailable")
print(f"  modality at peak: nmodes={df.nmodes[pk]}, bic={df.bic[pk]:.2f} -> "
      f"{'UNIMODAL (near-critical, like HUMAN)' if df.bic[pk]<0 else 'BIMODAL (switch-like, like MOUSE)'}")
frac_bi=np.mean([1 for i in range(nb) if (df.bic[i]>0)])
print(f"  fraction bimodal bins: {frac_bi:.2f}")
print("\nTRIAD: human=unimodal@peak(near-critical); mouse=bimodal@peak(switch); marmoset=?")
