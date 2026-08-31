# ============================================================================
# Figure S13  --  embryo-of-origin decomposition
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM transition reading and the ESRRB marker are independent of cESFW
#   - the key circularity control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   embryo-of-origin decomposition
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_embryo_decomp.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Embryo-of-origin decomposition of the human ICM order parameter (Methods).
Tests whether the ICM 'single unresolved population' is a pooling artefact:
between-embryo variance fraction (eta^2) vs a permutation null, and modality of
embryo-mean order parameters (unimodal continuum vs two clusters)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, sys, os
from sklearn.mixture import GaussianMixture
sys.path.insert(0, os.getcwd()); from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]; HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["a"]=ann[m]
ad.obs["emb"]=SI["Embryo"].astype(str).values[m]; ad.obs["ds"]=SI["Dataset"].astype(str).values[m]; ad.obs["day"]=SI["EmbryoDay"].astype(str).values[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="_e"); sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="_h")
phi1=z(al["NANOG"].values)-z(al["GATA6"].values); phi3=z(ad.obs["_e"].values)-z(ad.obs["_h"].values)
icm=ad.obs["a"].values=="ICM"; emb=ad.obs["emb"].values[icm]; ds=ad.obs["ds"].values[icm]; day=ad.obs["day"].values[icm]
def eta2(y,g):
    df=pd.DataFrame({"y":y,"g":g}); gm=df.y.mean()
    ssb=df.groupby("g").apply(lambda d:len(d)*(d.y.mean()-gm)**2).sum(); return ssb/((df.y-gm)**2).sum()
def dbic(x):
    X=np.asarray(x).reshape(-1,1); return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)
rng=np.random.RandomState(0)
print(f"ICM n={icm.sum()} from {len(set(emb))} embryos; median {int(pd.Series(emb).value_counts().median())} cells/embryo")
for lab,p in [("phi1 (NANOG-GATA6)",phi1[icm]),("phi3 (multigene)",phi3[icm])]:
    obs=eta2(p,emb); null=np.array([eta2(p,rng.permutation(emb)) for _ in range(1000)]); pv=(null>=obs).mean()
    em=pd.Series(p).groupby(emb).mean().values
    print(f"{lab}: eta2_embryo={obs:.2f} (null {null.mean():.2f}, p={pv:.3f}); eta2_dataset={eta2(p,ds):.2f}; eta2_day={eta2(p,day):.2f}; embryo-mean dBIC={dbic(em):+.1f}")
