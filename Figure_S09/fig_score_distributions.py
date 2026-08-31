# ============================================================================
# Figure S9  --  score distributions
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   ICM order-parameter distributions under the multigene epiblast-hypoblast
#   score (companion to Figure 6).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   score distributions
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S9. ICM order-parameter distributions under the multigene
#   epiblast-hypoblast score. Companion to Figure 6, using a multigene
#   epiblast-minus-hypoblast signature score in place of the two-gene
#   NANOG–GATA6 toggle. Under this composite the human ICM reads as a single
#   population (ΔBIC=-9.8), as does the marmoset (-9.6), while the mouse ICM
#   remains decisively bimodal (+138). That the human is decisively unimodal
#   under the multigene composite but only borderline under the two-gene
#   toggle (Figure 6), while the mouse is decisively bimodal under both,
#   indicates that the mixture interpretation of the human ICM is order-
#   parameter-dependent and not robustly supported. The mouse-versus-primate
#   contrast, by contrast, is robust to the choice of order parameter.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_score_distributions.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Supplementary companion to Fig 6: the ICM order-parameter distributions under
the MULTIGENE epiblast-minus-hypoblast SCORE (phi3), instead of the 2-gene
NANOG-GATA6 toggle. Shows the human ICM reads as a single population under the
multigene score, while the mouse ICM remains resolved (bimodal)."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]; HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
def bic(x): X=np.asarray(x).reshape(-1,1); return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)
def score_phi(ad):
    sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="_e")
    sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="_h")
    return z(ad.obs["_e"].values)-z(ad.obs["_h"].values)

# HUMAN
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["a"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
pH=score_phi(ad)[ad.obs["a"].values=="ICM"]
# MARMOSET
D="data/nhp/marmoset/RawCounts.csv"; meta=pd.read_csv(D,index_col=0,nrows=2); lin=meta.iloc[1]
cnt=pd.read_csv(D,skiprows=[1,2],index_col=0); cnt=cnt[~cnt.index.duplicated()]
cells=[c for c in cnt.columns if lin.get(c) in ("ICM_CS3","Epi_CS3","Hyp_CS3")]
adm=sc.AnnData(X=cnt[cells].T.values.astype(float)); adm.var_names=cnt.index; adm.obs_names=cells; adm.obs["l"]=[lin[c] for c in cells]
sc.pp.normalize_total(adm,target_sum=1e4); sc.pp.log1p(adm)
pM=score_phi(adm)[np.array(adm.obs["l"])=="ICM_CS3"]
# MOUSE
DATA="data/mouse_dev"; F={"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
                          "E4.5":["GSM3494336_Lib1-1_E4.5_counts.csv.gz","GSM3494337_Lib1-2_E4.5_counts.csv.gz"]}
fr,st=[],[]
for s,fs in F.items():
    for f in fs:
        d=pd.read_csv(os.path.join(DATA,f),index_col=0); d.index=[f"{f[:10]}_{i}" for i in range(len(d))]; fr.append(d); st+=[s]*len(d)
mt=pd.concat(fr).fillna(0); adk=sc.AnnData(X=mt.values.astype(float)); adk.var_names=[g.upper() for g in mt.columns]
adk.var_names_make_unique(); adk.obs_names=mt.index; adk.obs["stage"]=st
sc.pp.filter_cells(adk,min_genes=1500); sc.pp.normalize_total(adk,target_sum=1e4); sc.pp.log1p(adk)
sc.tl.score_genes(adk,[g for g in ["CDX2","GATA3","KRT8","KRT18","TEAD4","EOMES"] if g in adk.var_names],score_name="TE")
pK=score_phi(adk)[(adk.obs["TE"].values<0.15)&(adk.obs["stage"].values=="E3.5")]

data=[("Human ICM","#C0392B",pH),("Marmoset ICM","#E8873A",pM),("Mouse ICM (E3.5)","#2E86C1",pK)]
fig,ax=plt.subplots(1,3,figsize=(12,3.9),gridspec_kw={"wspace":0.3})
for i,(nm,c,x) in enumerate(data):
    b=bic(x); g=np.linspace(x.min()-.7,x.max()+.7,200); d=gaussian_kde(x)(g)
    ax[i].hist(x,bins=16,density=True,color=c,alpha=0.35,edgecolor="white",linewidth=0.5); ax[i].plot(g,d,color=c,lw=2.4)
    v="unimodal (single population)" if b<0 else "bimodal (resolved)"
    ax[i].set_title(f"{nm}  (n={len(x)})\n{v}  ($\\Delta$BIC={b:+.0f})",fontsize=10)
    ax[i].set_xlabel("$\\phi_3$ = multigene epiblast $-$ hypoblast score"); ax[i].set_ylabel("density" if i==0 else "")
    ax[i].axvline(0,color="#999",ls=":",lw=1); ax[i].spines[["top","right"]].set_visible(False)
fig.suptitle("ICM order-parameter distributions under the multigene epiblast$-$hypoblast score\n"
             "the human ICM reads as a single population under the multigene score; the mouse ICM remains resolved (bimodal)",fontsize=11,y=1.06)
out="plots/revision_results/Fig_score_distributions.png"; fig.savefig(out,dpi=200,bbox_inches="tight"); print("Saved:",out)
for nm,c,x in data: print(f"  {nm}: n={len(x)} bic={bic(x):+.2f}")