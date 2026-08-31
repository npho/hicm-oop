# ============================================================================
# Figure S10  --  landscape signatures under phi3 vs phi1
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM's unimodality is robust to the order parameter, but the
#   flattened/high-variance signature is NOT - an explicitly negative
#   robustness result.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   landscape signatures under phi3 vs phi1
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S10. The ICM’s unimodality is robust to the order parameter; the
#   flattened/high-variance signature is not. Robustness of the quasi-
#   potential signatures (Figure 5) to the choice of order parameter,
#   comparing the canonical two-gene toggle ϕ1=z(NANOG)-z(GATA6)with the
#   dropout-insensitive multigene composite ϕ3. (A) Variance of the order
#   parameter along diffusion pseudotime (normalised to each curve’s
#   maximum). Under ϕ1the variance peaks at the ICM; under ϕ3it continues to
#   rise through the ICM to the resolved hypoblast, so the variance maximum
#   is not localised to the ICM. (B) Dispersion of the ICM relative to the
#   resolved lineages: under ϕ1the ICM is ∼3-fold broader than the committed
#   epiblast/hypoblast (consistent with basin flattening), but under ϕ3this
#   reverses and the ICM is in fact tighter than the resolved lineages. The
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_landscape_robust.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Supplementary disclosure figure (Fig S12): robustness of the quasi-potential
signatures to the order parameter. The ICM's UNIMODALITY is robust; the
'flattened/high-variance' signature is specific to the two-gene toggle phi1
(inflated by NANOG dropout) and does not survive the multigene composite phi3."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]; HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]; ROOT="ICM/TE Branch"
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["a"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
sc.pp.neighbors(ad,n_neighbors=15,use_rep="X",metric="correlation"); sc.tl.diffmap(ad)
ad.uns["iroot"]=int(np.where(ad.obs["a"].values==ROOT)[0][0]); sc.tl.dpt(ad)
pt=ad.obs["dpt_pseudotime"].values; a=ad.obs["a"].values
sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="_e"); sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="_h")
phi1=z(al["NANOG"].values)-z(al["GATA6"].values); phi3=z(ad.obs["_e"].values)-z(ad.obs["_h"].values)
nb=8; edges=np.quantile(pt,np.linspace(0,1,nb+1)); edges[-1]+=1e-9
mid=[pt[(pt>=edges[i])&(pt<edges[i+1])].mean() for i in range(nb)]
stg=[pd.Series(a[(pt>=edges[i])&(pt<edges[i+1])]).mode().iloc[0] for i in range(nb)]
def vartraj(phi): return [np.var(phi[(pt>=edges[i])&(pt<edges[i+1])]) for i in range(nb)]
v1=vartraj(phi1); v3=vartraj(phi3)

fig,ax=plt.subplots(1,2,figsize=(12,4.4),gridspec_kw={"wspace":0.32,"width_ratios":[1.35,1]})
# Panel A: variance trajectory (normalised to own max) under both order parameters
axA=ax[0]
axA.plot(mid,np.array(v1)/max(v1),"o-",color="#C0392B",lw=2.4,ms=6,label="$\\phi_1$ = z(NANOG)$-$z(GATA6)  (2-gene)")
axA.plot(mid,np.array(v3)/max(v3),"s--",color="#6C3483",lw=2.4,ms=6,label="$\\phi_3$ = multigene epi$-$hypo")
icm_bin=[i for i in range(nb) if stg[i]=="ICM"][0]
axA.axvline(mid[icm_bin],color="#888",ls=":",lw=1); axA.text(mid[icm_bin],1.03,"ICM",ha="center",fontsize=8.5,color="#555")
axA.set_xlabel("diffusion pseudotime"); axA.set_ylabel("Var(order parameter), normalised")
axA.set_title("A  Variance peaks at the ICM only under the 2-gene toggle\n($\\phi_3$ keeps rising to the resolved hypoblast)",fontsize=9.8)
axA.legend(fontsize=8,frameon=False,loc="lower right"); axA.spines[["top","right"]].set_visible(False)
# Panel B: ICM vs committed dispersion under each order parameter
axB=ax[1]
comm1=np.mean([np.var(phi1[a=="preIm-Epi"]),np.var(phi1[a=="Hypo"])])
comm3=np.mean([np.var(phi3[a=="preIm-Epi"]),np.var(phi3[a=="Hypo"])])
icm1=np.var(phi1[a=="ICM"]); icm3=np.var(phi3[a=="ICM"])
xp=np.arange(2); w=0.36
axB.bar(xp-w/2,[icm1,icm3],w,color="#E59866",label="ICM")
axB.bar(xp+w/2,[comm1,comm3],w,color="#85929E",label="resolved lineages (ref)")
for i,(a1,b1) in enumerate([(icm1,comm1),(icm3,comm3)]):
    axB.text(i,max(a1,b1)+0.08,("ICM broader" if a1>b1 else "ICM tighter"),ha="center",fontsize=8.5,
             color=("#C0392B" if a1>b1 else "#6C3483"),fontweight="bold")
axB.set_xticks(xp); axB.set_xticklabels(["$\\phi_1$ (2-gene)","$\\phi_3$ (multigene)"])
axB.set_ylabel("Var of order parameter"); axB.legend(fontsize=8,frameon=False)
axB.set_title("B  'ICM more dispersed than resolved' is $\\phi_1$-specific\n(reverses under $\\phi_3$)",fontsize=9.8)
axB.spines[["top","right"]].set_visible(False)
fig.suptitle("The ICM's unimodality is robust to the order parameter; the flattened/high-variance signature is not\n"
             "(NANOG is undetected in ~16% of ICM cells, inflating the two-gene $\\phi_1$ variance at the ICM)",fontsize=11,y=1.04)
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"plots","revision_results")
out=os.path.join(OUT,"Fig_landscape_robust.png"); fig.savefig(out,dpi=200,bbox_inches="tight",facecolor="white")
print("Saved:",out); print("phi1: ICM var %.2f vs committed %.2f"%(icm1,comm1)); print("phi3: ICM var %.2f vs committed %.2f"%(icm3,comm3))
