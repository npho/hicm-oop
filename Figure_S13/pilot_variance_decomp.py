# ============================================================================
# Figure S13  --  variance decomposition vs technical covariates
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM transition reading and the ESRRB marker are independent of cESFW
#   - the key circularity control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   variance decomposition vs technical covariates
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/pilot_variance_decomp.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""HONEST pilot: how much of the human ICM order-parameter (phi) heterogeneity is
TECHNICAL (cell-cycle, sequencing depth, dataset-of-origin) vs biological?
Reported whatever it shows. NOT a 'purify until near-critical wins' analysis."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, statsmodels.formula.api as smf
from sklearn.mixture import GaussianMixture
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
def dbic(x):
    X=np.asarray(x).reshape(-1,1)
    return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)

LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
raw=C.values[m].astype(float)
ad=sc.AnnData(X=raw.copy()); ad.var_names=C.columns; ad.obs_names=C.index[m]
ad.obs["ann"]=ann[m]; ad.obs["ds"]=SI["Dataset"].astype(str).values[m]
# depth from raw counts BEFORE normalisation
ad.obs["logdepth"]=np.log1p(raw.sum(1)); ad.obs["ngenes"]=(raw>0).sum(1)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
# cell-cycle score (Tirosh subset present in data)
Sg=[g for g in ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"] if g in ad.var_names]
Gg=[g for g in ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"] if g in ad.var_names]
sc.tl.score_genes(ad,Sg,score_name="Sscore"); sc.tl.score_genes(ad,Gg,score_name="G2Mscore")

phi=z(al["NANOG"].values)-z(al["GATA6"].values)
ad.obs["phi"]=phi
icm=ad.obs[ad.obs["ann"].values=="ICM"].copy()
print("human ICM cells: %d | var(phi)=%.3f | DeltaBIC=%.2f (%s)"%(
      len(icm),icm["phi"].var(),dbic(icm["phi"].values),"unimodal" if dbic(icm["phi"].values)<0 else "bimodal"))

print("\n=== variance of ICM phi explained (R^2) by each confounder alone ===")
for f in ["Sscore","G2Mscore","logdepth","ngenes","C(ds)"]:
    try:
        r=smf.ols("phi ~ "+f,data=icm).fit(); print("  %-12s R^2 = %.3f"%(f,r.rsquared))
    except Exception as e: print("  %-12s failed (%s)"%(f,e))
full=smf.ols("phi ~ Sscore + G2Mscore + logdepth + C(ds)",data=icm).fit()
print("  %-12s R^2 = %.3f  (all technical covariates jointly)"%("FULL",full.rsquared))
print("\n=> %.0f%% of ICM phi variance is attributable to the measured technical covariates; %.0f%% is residual/biological."
      %(100*full.rsquared,100*(1-full.rsquared)))
print("\n=== diagnostic only (NOT the headline): modality of residual phi after removing technical variance ===")
res=full.resid.values
print("  residual phi: var=%.3f  DeltaBIC=%.2f -> %s"%(np.var(res),dbic(res),"unimodal" if dbic(res)<0 else "bimodal"))
print("  (report honestly whichever way it goes; do not use to overturn the ambiguity)")