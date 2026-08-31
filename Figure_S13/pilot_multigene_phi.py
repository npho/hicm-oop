# ============================================================================
# Figure S13  --  multigene order parameter per dataset
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM transition reading and the ESRRB marker are independent of cESFW
#   - the key circularity control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   multigene order parameter per dataset
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S13. The ICM transition reading and the ESRRB ICM-stage marker
#   are both independent of cESFW. Modality and marker analyses repeated
#   with cESFW removed from every step (standard normalisation, canonical
#   markers, full transcriptome; cells defined by independent Stirparo et
#   al. 2018 labels). (A) φ = z(epiblast) − z(hypoblast): the ICM (∪ early
#   ICM) is a single, unimodal population positioned between the epiblast
#   and hypoblast poles. (B) Positive control (n matched): a synthetic 50/50
#   mixture of the independently annotated epiblast and hypoblast cells is
#   multimodal (dip p<0.001) whereas the ICM is unimodal (dip p=0.96); a
#   stage-only gate (embryonic day 5 inner cells, no lineage presupposition)
#   gives the same result (dip p=0.98). (C) ESRRB, high in the ICM/early
#   compartment, is extinguished in the naïve epiblast while NANOG and KLF17
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/pilot_multigene_phi.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Does the human ICM classification change under a MULTIGENE order parameter
(vs the 2-gene NANOG-GATA6)? Redo the two checks that drove the deflation, the
BOOTSTRAP three-model win-frequency and the PER-DATASET / annotation modality,
on phi3 = multigene(epiblast) - multigene(hypoblast). Report honestly."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]; HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
def bic_gmm(x,k): X=np.asarray(x).reshape(-1,1); return GaussianMixture(k,n_init=3,random_state=0).fit(X).bic(X)
def bic_att(x,vr): x=np.asarray(x); return -2*norm.logpdf(x,x.mean(),np.sqrt(vr)).sum()+np.log(len(x))
def dbic(x): return bic_gmm(x,1)-bic_gmm(x,2)
def boot(icm,vr,B=500):
    w=Counter(); rng=np.random.RandomState(0)
    for _ in range(B):
        x=np.asarray(icm)[rng.randint(0,len(icm),len(icm))]
        b=[("attr",bic_att(x,vr)),("mix",bic_gmm(x,2)),("near",bic_gmm(x,1))]
        w[min(b,key=lambda t:t[1])[0]]+=1
    return {k:round(w.get(k,0)/B,3) for k in ["attr","mix","near"]}

LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]
ad.obs["a"]=ann[m]; ad.obs["ds"]=SI["Dataset"].astype(str).values[m]; ad.obs["stir"]=SI["Stirparo_Labels"].astype(str).values[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="ep")
sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="hy")
phi1=z(al["NANOG"].values)-z(al["GATA6"].values)
phi3=z(ad.obs["ep"].values)-z(ad.obs["hy"].values)
aH=ad.obs["a"].values

for nm,phi in [("phi1 = NANOG-GATA6 (2-gene)",phi1),("phi3 = multigene epi-hyp score",phi3)]:
    icm=phi[aH=="ICM"]; vr=float(np.mean([np.var(phi[aH=="preIm-Epi"]),np.var(phi[aH=="Hypo"])]))
    print("\n===== %s ====="%nm)
    print("  ICM: n=%d var=%.2f  DeltaBIC=%.2f (%s)"%(len(icm),np.var(icm),dbic(icm),"unimodal" if dbic(icm)<0 else "bimodal"))
    print("  bootstrap three-model win-freq:", boot(icm,vr))
    print("  per-dataset DeltaBIC:")
    for d in sorted(set(ad.obs["ds"].values)):
        x=phi[(aH=="ICM")&(ad.obs["ds"].values==d)]
        if len(x)>=12: print("     %-16s n=%3d  DeltaBIC=%6.2f -> %s"%(d,len(x),dbic(x),"unimodal" if dbic(x)<0 else "bimodal"))
    xs=phi[np.isin(ad.obs["stir"].values,["ICM","Early ICM"])]
    print("  Stirparo ICM (n=%d): DeltaBIC=%.2f -> %s"%(len(xs),dbic(xs),"unimodal" if dbic(xs)<0 else "bimodal"))