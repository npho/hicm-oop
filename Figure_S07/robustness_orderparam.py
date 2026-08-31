# ============================================================================
# Figure S7  --  three order parameters compared
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mouse-versus-primate ICM contrast is robust to the choice of order
#   parameter - guards against the result being an artefact of NANOG/GATA6.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   three order parameters compared
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S7. The mouse-versus-primate ICM contrast is robust to the order
#   parameter; the strength of the human classification is not. Two-
#   component Gaussian-mixture support (ΔBIC; >0bimodal, <0unimodal) for the
#   ICM of each species, computed with three order parameters:
#   ϕ1=z(NANOG)-z(GATA6)(the two-gene toggle used in the main text), ϕ2=the
#   first principal component of the lineage genes (epiblast and hypoblast
#   markers, standardised and sign-aligned), and ϕ3=a multigene epiblast-
#   minus-hypoblast score. The mouse ICM is strongly bimodal (+87 to +138)
#   and both primate ICMs unimodal (ΔBIC <0) at the point estimate under all
#   three, so the cross-species contrast does not rest on any single marker
#   pair. The human point estimate, however, is only weakly unimodal under
#   the canonical two-gene toggle (-2.7), though stronger under the
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/robustness_orderparam.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""ROBUSTNESS: does the near-critical (unimodal primate) vs switch-like (bimodal
mouse) ICM result survive alternative definitions of the order parameter phi?
phi1 = z(NANOG) - z(GATA6)  [used in the paper]
phi2 = PC1 of lineage genes (epiblast + hypoblast markers), sign-aligned
phi3 = multigene epiblast score - hypoblast score (scanpy score_genes)
Test: DeltaBIC (1- vs 2-component GMM) of the ICM population, per species per phi."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z = lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=["NANOG","SOX2","POU5F1","KLF17","TDGF1","DPPA5","PRDM14"]
HYP=["GATA6","GATA4","SOX17","PDGFRA","FOXA2","HNF1B"]
def bic(x):
    X=np.asarray(x).reshape(-1,1)
    return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)

def three_phi(al, ad, icm_mask):
    """compute phi1/phi2/phi3 over all compartment cells; return ICM-subset values."""
    epi=[g for g in EPI if g in al.columns]; hyp=[g for g in HYP if g in al.columns]
    # phi1
    phi1 = z(al["NANOG"].values) - z(al["GATA6"].values)
    # phi3: multigene score difference
    sc.tl.score_genes(ad, epi, score_name="_EPI", ctrl_size=50)
    sc.tl.score_genes(ad, hyp, score_name="_HYP", ctrl_size=50)
    phi3 = z(ad.obs["_EPI"].values) - z(ad.obs["_HYP"].values)
    # phi2: PC1 of lineage genes (standardised)
    genes=epi+hyp; M=al[genes].values
    M=(M-M.mean(0))/(M.std(0)+1e-9)
    pc1=PCA(n_components=1, random_state=0).fit_transform(M)[:,0]
    if np.corrcoef(pc1, phi1)[0,1] < 0: pc1=-pc1     # sign-align to epiblast-positive
    phi2=z(pc1)
    return {"phi1 = z(NANOG)-z(GATA6)":phi1[icm_mask],
            "phi2 = PC1(lineage genes)":phi2[icm_mask],
            "phi3 = epi-score - hyp-score":phi3[icm_mask]}

# ---------- HUMAN ----------
LIN=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
SI,C,cG,E,*_=load_embryo_data(); ann=SI["Manual_Annotations"].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs["a"]=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
H=three_phi(al, ad, ad.obs["a"].values=="ICM")

# ---------- MARMOSET ----------
D="data/nhp/marmoset/RawCounts.csv"; meta=pd.read_csv(D,index_col=0,nrows=2); lineage=meta.iloc[1]
counts=pd.read_csv(D,skiprows=[1,2],index_col=0); counts=counts[~counts.index.duplicated()]
cells=[c for c in counts.columns if lineage.get(c) in ("ICM_CS3","Epi_CS3","Hyp_CS3")]
adm=sc.AnnData(X=counts[cells].T.values.astype(float)); adm.var_names=counts.index; adm.obs_names=cells
adm.obs["lin"]=[lineage[c] for c in cells]
sc.pp.normalize_total(adm,target_sum=1e4); sc.pp.log1p(adm); alm=pd.DataFrame(adm.X,index=adm.obs_names,columns=adm.var_names)
Ma=three_phi(alm, adm, np.array(adm.obs["lin"])=="ICM_CS3")

# ---------- MOUSE ----------
DATA="data/mouse_dev"
F={"E3.5":["GSM3494334_Lib1-3_E3.5_counts.csv.gz","GSM3494335_Lib1-4_E3.5_counts.csv.gz"],
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
Mo=three_phi(alk, adk, inner)

species=[("Human",H,"#C0392B"),("Marmoset",Ma,"#E8873A"),("Mouse",Mo,"#2E86C1")]
phis=list(H.keys())
print("\n=== ROBUSTNESS: DeltaBIC of ICM (>0 bimodal/switch, <0 unimodal/near-critical) ===")
print("%-14s %14s %14s %14s"%("species","phi1(2-gene)","phi2(PC1)","phi3(scores)"))
for nm,dct,_ in species:
    print("%-14s %14.2f %14.2f %14.2f"%(nm, bic(dct[phis[0]]), bic(dct[phis[1]]), bic(dct[phis[2]])))

# figure: grouped bars
fig,ax=plt.subplots(figsize=(8.2,4.4)); x=np.arange(3); w=0.26
for j,(nm,dct,c) in enumerate(species):
    vals=[bic(dct[p]) for p in phis]
    ax.bar(x+(j-1)*w, vals, w, color=c, alpha=0.9, edgecolor="white", label=nm)
ax.axhline(0,color="#333",lw=1)
ax.set_xticks(x); ax.set_xticklabels(["$\\phi_1$: z(NANOG)-z(GATA6)","$\\phi_2$: PC1(lineage)","$\\phi_3$: epi-hyp score"],fontsize=9)
ax.set_ylabel("$\\Delta$BIC of ICM  ( >0 bimodal / <0 unimodal )")
ax.set_title("The near-critical (unimodal) primate ICM vs switch-like (bimodal) mouse ICM\nis robust to the definition of the order parameter",fontsize=11)
ax.legend(frameon=False,fontsize=9); ax.spines[["top","right"]].set_visible(False)
ax.text(0.01,0.02,"primates below 0 (unimodal) under all three definitions; mouse above 0 (bimodal)",
        transform=ax.transAxes,fontsize=8,color="#666")
OUT="plots/revision_results"; os.makedirs(OUT,exist_ok=True)
out=os.path.join(OUT,"Fig_robustness_orderparam.png"); fig.savefig(out,dpi=200,bbox_inches="tight")
print("\nSaved:",out)
