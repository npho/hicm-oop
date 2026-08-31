# ============================================================================
# Figure 6  --  A: modality of the ICM order parameter, both order parameters
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mouse ICM is already resolved into two populations; the human and
#   marmoset ICM are not. The species contrast that anchors the transition
#   claim.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A: modality of the ICM order parameter, both order parameters
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/theory_crossspecies_both.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Merged cross-species Figure 6 (v46): the ICM order-parameter distribution for
three species under BOTH order parameters, with the dropout-robust MULTIGENE
epiblast-minus-hypoblast composite (phi3) as the primary row and the canonical
two-gene NANOG-GATA6 toggle (phi1) as validation.

Message: under the robust composite the human and marmoset ICMs are decisively a
single unresolved population while the mouse ICM is already resolved into two;
the human two-gene toggle looks borderline only because NANOG is undetected in a
sizeable fraction of ICM cells, exactly as expected of an unresolved transition."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z = lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]
HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
def bic(x):
    X=np.asarray(x).reshape(-1,1)
    return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)
def score_phi(ad):
    sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name="_e")
    sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name="_h")
    return z(ad.obs["_e"].values)-z(ad.obs["_h"].values)

# ---------------- HUMAN (cESFW atlas 'ICM') ----------------
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["a"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
icm=ad.obs["a"].values=="ICM"
phi1_H=(z(al["NANOG"].values)-z(al["GATA6"].values))[icm]
phi3_H=score_phi(ad)[icm]
# NANOG dropout fraction within ICM (raw counts)
nanog_raw=C.values[m][:, list(C.columns).index("NANOG")][icm]
gata6_raw=C.values[m][:, list(C.columns).index("GATA6")][icm]
drop_nanog=float(np.mean(nanog_raw==0)); drop_gata6=float(np.mean(gata6_raw==0))

# ---------------- MARMOSET (Boroviak ICM_CS3) ----------------
D="data/nhp/marmoset/RawCounts.csv"; meta=pd.read_csv(D,index_col=0,nrows=2); lin=meta.iloc[1]
cnt=pd.read_csv(D,skiprows=[1,2],index_col=0); cnt=cnt[~cnt.index.duplicated()]
cells=[c for c in cnt.columns if lin.get(c) in ("ICM_CS3","Epi_CS3","Hyp_CS3")]
adm=sc.AnnData(X=cnt[cells].T.values.astype(float)); adm.var_names=cnt.index; adm.obs_names=cells; adm.obs["l"]=[lin[c] for c in cells]
sc.pp.normalize_total(adm,target_sum=1e4); sc.pp.log1p(adm)
alm=pd.DataFrame(adm.X,index=adm.obs_names,columns=adm.var_names)
mi=np.array(adm.obs["l"])=="ICM_CS3"
phi1_M=(z(alm["NANOG"].values)-z(alm["GATA6"].values))[mi]
phi3_M=score_phi(adm)[mi]

# ---------------- MOUSE (Nowotschin E3.5 inner) ----------------
DATA="data/mouse_dev"; F={"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
                          "E4.5":["GSM3494336_Lib1-1_E4.5_counts.csv.gz","GSM3494337_Lib1-2_E4.5_counts.csv.gz"]}
fr,st=[],[]
for s,fs in F.items():
    for f in fs:
        d=pd.read_csv(os.path.join(DATA,f),index_col=0); d.index=[f"{f[:10]}_{i}" for i in range(len(d))]; fr.append(d); st+=[s]*len(d)
mt=pd.concat(fr).fillna(0); adk=sc.AnnData(X=mt.values.astype(float)); adk.var_names=[g.upper() for g in mt.columns]
adk.var_names_make_unique(); adk.obs_names=mt.index; adk.obs["stage"]=st
sc.pp.filter_cells(adk,min_genes=1500); sc.pp.normalize_total(adk,target_sum=1e4); sc.pp.log1p(adk)
alk=pd.DataFrame(adk.X,index=adk.obs_names,columns=adk.var_names)
sc.tl.score_genes(adk,[g for g in ["CDX2","GATA3","KRT8","KRT18","TEAD4","EOMES"] if g in adk.var_names],score_name="TE")
inner=(adk.obs["TE"].values<0.15)&(adk.obs["stage"].values=="E3.5")
phi1_K=(z(alk["NANOG"].values)-z(alk["GATA6"].values))[inner]
phi3_K=score_phi(adk)[inner]

SPEC=[("Human","#C0392B",phi3_H,phi1_H),("Marmoset","#E8873A",phi3_M,phi1_M),("Mouse (E3.5)","#2E86C1",phi3_K,phi1_K)]
OUT=os.path.join("plots","revision_results"); os.makedirs(OUT,exist_ok=True)

def verdict(b):
    return "unimodal (single population)" if b<-5 else "bimodal (resolved / switch-like)" if b>5 else "borderline (unresolved)"

fig=plt.figure(figsize=(15,7.6))
gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,0.95],height_ratios=[1,1],wspace=0.30,hspace=0.52)
rows=[("$\\phi_3$ = multigene epiblast $-$ hypoblast score",2,"ROBUST (dropout-insensitive)"),
      ("$\\phi_1$ = z(NANOG) $-$ z(GATA6)",3,"canonical two-gene toggle")]
b3=[]; b1=[]
for r,(xlab,idx,tag) in enumerate(rows):
    for cidx,(nm,c,x3,x1) in enumerate(SPEC):
        x=x3 if r==0 else x1
        ax=fig.add_subplot(gs[r,cidx])
        b=bic(x); (b3 if r==0 else b1).append(b)
        g=np.linspace(x.min()-.7,x.max()+.7,200); d=gaussian_kde(x)(g)
        ax.hist(x,bins=16,density=True,color=c,alpha=0.32,edgecolor="white",linewidth=0.5)
        ax.plot(g,d,color=c,lw=2.4)
        ax.set_title(f"{nm} ICM  (n={len(x)})\n{verdict(b)}   ($\\Delta$BIC={b:+.0f})",fontsize=9.6)
        ax.set_xlabel(xlab,fontsize=9); ax.set_ylabel("density" if cidx==0 else "")
        ax.axvline(0,color="#999",ls=":",lw=1); ax.spines[["top","right"]].set_visible(False)
        if cidx==0:
            ax.text(-0.30,0.5,("PRIMARY\n"+tag) if r==0 else ("VALIDATION\n"+tag),transform=ax.transAxes,
                    rotation=90,ha="center",va="center",fontsize=8.6,color="#444",fontweight="bold")

# ---- summary panel (spans both rows, right column) ----
axs=fig.add_subplot(gs[:,3])
names=[s[0].split(" (")[0] for s in SPEC]; xpos=np.arange(len(names)); w=0.38
from matplotlib.patches import Patch
axs.bar(xpos-w/2,b3,w,color=[s[1] for s in SPEC],alpha=0.95,edgecolor="white")
axs.bar(xpos+w/2,b1,w,color=[s[1] for s in SPEC],alpha=0.45,edgecolor="white",hatch="//")
axs.axhline(0,color="#333",lw=1); axs.axhspan(-5,5,color="#ccc",alpha=0.35,zorder=0)
axs.set_xticks(xpos); axs.set_xticklabels(names,rotation=18,fontsize=9)
axs.set_ylabel("$\\Delta$BIC   ( >0 bimodal / resolved )",fontsize=9.5)
axs.set_title("ICM modality across species\n(both order parameters)",fontsize=10)
axs.legend(handles=[Patch(facecolor="#777",edgecolor="white",label="$\\phi_3$ multigene (robust)"),
                    Patch(facecolor="#bbb",edgecolor="white",hatch="//",label="$\\phi_1$ 2-gene toggle")],
           fontsize=8,frameon=False,loc="upper left")
axs.spines[["top","right"]].set_visible(False)
axs.text(0.5,-0.10,"grey band = inconclusive zone.\nMouse resolved under both; primates single\npopulation (decisively under $\\phi_3$).",
         transform=axs.transAxes,ha="center",va="top",fontsize=7.6,color="#666")

fig.suptitle("The mouse ICM is resolved into two populations; the primate ICM is a single, unresolved population\n"
             "Robust under the dropout-insensitive multigene composite ($\\phi_3$, top); the human two-gene toggle ($\\phi_1$, bottom) looks borderline because NANOG is undetected in %.0f%% of human ICM cells"%(100*drop_nanog),
             fontsize=12,y=0.99)
out=os.path.join(OUT,"Fig_crossspecies_both.png"); fig.savefig(out,dpi=200,bbox_inches="tight",facecolor="white")
print("Saved:",out)
print("=== NUMBERS (for text consistency) ===")
for r,(xlab,idx,tag) in enumerate(rows):
    lab="phi3 (multigene)" if r==0 else "phi1 (2-gene)"
    for cidx,(nm,c,x3,x1) in enumerate(SPEC):
        x=x3 if r==0 else x1
        print(f"  {lab:20s} {nm:14s} n={len(x):4d}  var={np.var(x):5.2f}  dBIC={bic(x):+7.2f}")
print(f"  NANOG dropout in human ICM: {100*drop_nanog:.1f}%   GATA6 dropout: {100*drop_gata6:.1f}%")
