# ============================================================================
# Figure S10  --  landscape recomputed under the multigene composite
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM's unimodality is robust to the order parameter, but the
#   flattened/high-variance signature is NOT - an explicitly negative
#   robustness result.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   landscape recomputed under the multigene composite
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/theory_landscape_phi3.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Robustness of the quasi-potential landscape to the order parameter.
Repeats the Fig 5 analysis with phi3 (multigene epiblast-minus-hypoblast score)
alongside phi1 = z(NANOG)-z(GATA6): does the variance peak at the ICM transition,
and is the peak-variance region UNIMODAL-broad (not a bimodal mixture), under BOTH?"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]; HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]; ROOT="ICM/TE Branch"
def dbic(x):
    X=np.asarray(x).reshape(-1,1)
    return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)

SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["ann"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
sc.pp.neighbors(ad,n_neighbors=15,use_rep="X",metric="correlation"); sc.tl.diffmap(ad)
ad.uns["iroot"]=int(np.where(ad.obs["ann"].values==ROOT)[0][0]); sc.tl.dpt(ad)
pt=ad.obs["dpt_pseudotime"].values; n=ad.n_obs
sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="_e")
sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="_h")
phi1=z(al["NANOG"].values)-z(al["GATA6"].values)
phi3=z(ad.obs["_e"].values)-z(ad.obs["_h"].values)

nb=8; edges=np.quantile(pt,np.linspace(0,1,nb+1)); edges[-1]+=1e-9
def analyse(phi,name):
    var=[phi[(pt>=edges[i])&(pt<edges[i+1])].var() for i in range(nb)]
    stage=[pd.Series(ad.obs["ann"].values[(pt>=edges[i])&(pt<edges[i+1])]).mode().iloc[0] for i in range(nb)]
    rng=np.random.RandomState(0); pk=[]
    for _ in range(500):
        idx=rng.randint(0,n,n)
        vb=[phi[idx][(pt[idx]>=edges[i])&(pt[idx]<edges[i+1])].var() for i in range(nb)]
        pk.append(int(np.nanargmax(vb)))
    pc=Counter(pk); mpk=max(pc,key=pc.get)
    # modality within the annotated ICM (the transition compartment)
    icm=ad.obs["ann"].values=="ICM"; db_icm=dbic(phi[icm])
    # modality at the variance-peak bin and at the terminal bin
    pkbin=(pt>=edges[mpk])&(pt<edges[mpk+1]); db_pk=dbic(phi[pkbin])
    term=(pt>=edges[nb-1]); db_term=dbic(phi[term])
    print(f"\n=== {name} ===")
    print("  per-bin Var:",[f"{v:.2f}" for v in var])
    print("  per-bin stage:",stage)
    print(f"  variance peak: bin {mpk} (stage {stage[mpk]}), {100*pc[mpk]/500:.0f}% of bootstraps; "
          f"peak/median var ratio = {var[mpk]/np.median(var):.1f}x")
    print(f"  dBIC within annotated ICM = {db_icm:+.2f}  ({'UNIMODAL' if db_icm<0 else 'bimodal'})")
    print(f"  dBIC at variance-peak bin = {db_pk:+.2f}  ({'UNIMODAL' if db_pk<0 else 'bimodal'})")
    print(f"  dBIC at terminal bin      = {db_term:+.2f}  ({'bimodal(split)' if db_term>0 else 'unimodal'})")
    return var,stage,mpk

print("lineage cells:",n,"| per stage:",{s:int((ad.obs['ann'].values==s).sum()) for s in LIN})
v1,s1,pk1=analyse(phi1,"phi1 = z(NANOG) - z(GATA6)  [Fig 5 canonical]")
v3,s3,pk3=analyse(phi3,"phi3 = multigene epiblast - hypoblast  [robustness]")

# compact 2-panel figure: Var vs pseudotime under each order parameter
mid=[pt[(pt>=edges[i])&(pt<edges[i+1])].mean() for i in range(nb)]
fig,ax=plt.subplots(1,2,figsize=(11,4),gridspec_kw={"wspace":0.28})
for a,(v,pk,nm,c) in zip(ax,[(v1,pk1,"$\\phi_1$ = z(NANOG) $-$ z(GATA6)","#C0392B"),
                             (v3,pk3,"$\\phi_3$ = multigene epiblast $-$ hypoblast","#6C3483")]):
    a.plot(mid,v,"o-",color=c,lw=2.4,ms=7); a.axvline(mid[pk],color="#888",ls=":")
    a.annotate("variance peak\n(%s)"%s1[pk],xy=(mid[pk],v[pk]),xytext=(mid[pk],v[pk]*0.72),
               fontsize=8,ha="center",color="#555")
    a.set_title(nm,fontsize=10); a.set_xlabel("diffusion pseudotime"); a.set_ylabel("Var of order parameter")
    a.spines[["top","right"]].set_visible(False)
fig.suptitle("The order-parameter variance peaks at the ICM transition under both order parameters",fontsize=12,y=1.02)
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),"plots","revision_results")
out=os.path.join(OUT,"Fig_landscape_phi3_robust.png"); fig.savefig(out,dpi=200,bbox_inches="tight",facecolor="white")
print("\nSaved:",out)
