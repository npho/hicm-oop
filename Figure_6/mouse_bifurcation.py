# ============================================================================
# Figure 6  --  mouse ICM bifurcation
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mouse ICM is already resolved into two populations; the human and
#   marmoset ICM are not. The species contrast that anchors the transition
#   claim.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   mouse ICM bifurcation
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/mouse_bifurcation.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""MOUSE critical-signature test (cross-species contrast with human).
Nowotschin 2019 E3.5->E4.5 (10x). phi = z(NANOG)-z(GATA6) on inner cells
(ICM/epi/PrE; TE excluded). KEY question: is mouse ALREADY BIMODAL at its
variance peak (switch-like / salt-and-pepper) vs human's unimodal-broad
near-critical window (lingering)?"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
np.seterr(all="ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE,"data","mouse_dev")

FILES = {"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
         "E4.5":["GSM3494336_Lib1-1_E4.5_counts.csv.gz","GSM3494337_Lib1-2_E4.5_counts.csv.gz"]}
frames, stages = [], []
for stg, fs in FILES.items():
    for f in fs:
        d = pd.read_csv(os.path.join(DATA,f), index_col=0); d.index=[f"{f[:10]}_{i}" for i in range(len(d))]
        frames.append(d); stages += [stg]*len(d)
mat = pd.concat(frames, axis=0).fillna(0)
ad = sc.AnnData(X=mat.values.astype(float)); ad.var_names=[g.upper() for g in mat.columns]
ad.var_names_make_unique(); ad.obs_names=mat.index; ad.obs["stage"]=stages
sc.pp.filter_cells(ad, min_genes=1500); sc.pp.filter_genes(ad, min_cells=3)
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
al = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)
def sc_score(genes, nm):
    g=[x for x in genes if x in ad.var_names]; sc.tl.score_genes(ad,g,score_name=nm,ctrl_size=50)
sc_score(["CDX2","GATA3","KRT8","KRT18","TEAD4","EOMES","KRT7"], "TE")
sc_score(["NANOG","SOX2","POU5F1","FGF4","ESRRB","ZFP42"], "EPI")
sc_score(["GATA6","GATA4","SOX17","PDGFRA","FOXA2","DAB2"], "PRE")
# inner cells = non-TE (ICM/epi/PrE); exclude trophectoderm
inner = ad.obs["TE"].values < 0.15
adi = ad[inner].copy(); ali = al.loc[adi.obs_names]
print(f"mouse inner cells: {adi.n_obs} (E3.5={int((adi.obs['stage']=='E3.5').sum())}, E4.5={int((adi.obs['stage']=='E4.5').sum())})")

z = lambda v:(v-v.mean())/(v.std()+1e-9)
phi = z(ali["NANOG"].values) - z(ali["GATA6"].values)
sc.pp.neighbors(adi, n_neighbors=15, use_rep="X", metric="correlation"); sc.tl.diffmap(adi)
adi.uns["iroot"]=int(np.where(adi.obs["stage"].values=="E3.5")[0][np.argmax(
    (z(ali["NANOG"].values)+z(ali["GATA6"].values))[adi.obs["stage"].values=="E3.5"])])  # root: co-expressing E3.5 cell
sc.tl.dpt(adi); pt = adi.obs["dpt_pseudotime"].values

def fita(x):
    if len(x)<20: return np.nan
    g=np.linspace(np.percentile(x,2)-.3,np.percentile(x,98)+.3,120); d=gaussian_kde(x)(g); U=-np.log(d+1e-6)
    w=d/d.sum(); Xb=np.vstack([np.ones_like(g),g,g**2,g**4]).T
    beta=np.linalg.lstsq((Xb*w[:,None]),U*w,rcond=None)[0]; return 2*beta[2]
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
    stg=pd.Series(adi.obs["stage"].values[sel]).mode()
    rows.append(dict(mid=pt[sel].mean(),nn=int(sel.sum()),stage=stg.iloc[0] if len(stg) else "?",
                     var=x.var(),a=fita(x),nmodes=nmodes(x,grid),bic=bic(x)))
df=pd.DataFrame(rows)
pk=int(np.nanargmax(df["var"].values))
print("\n=== MOUSE per pseudotime bin ==="); print(df.round(3).to_string(index=False))
print(f"\nvariance peak: bin {pk} (stage {df.stage[pk]}, pt={df.mid[pk]:.3f})")
print(f"  at peak: nmodes={df.nmodes[pk]}, bic={df.bic[pk]:.2f} -> "
      f"{'UNIMODAL (lingering/near-critical, like human)' if df.bic[pk]<0 else 'BIMODAL (switch-like / salt-and-pepper)'}")
firstbi=next((i for i in range(nb) if df.nmodes[i]>=2 and df.bic[i]>0), None)
print(f"  first BIMODAL bin: {firstbi} (fraction of trajectory: {df.mid[firstbi]:.2f})" if firstbi is not None else "  never clearly bimodal")
# overall: fraction of trajectory that is bimodal (switch-like index)
frac_bi = np.mean([1 for i in range(nb) if df.bic[i]>0])
print(f"  fraction of bins that are bimodal: {frac_bi:.2f}")
print(f"\nCONTRAST vs HUMAN (human: variance peak UNIMODAL, bic=-8.7; first bimodal at 75% of trajectory)")
