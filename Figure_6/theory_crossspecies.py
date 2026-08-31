# ============================================================================
# Figure 6  --  cross-species modality (single order parameter)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The mouse ICM is already resolved into two populations; the human and
#   marmoset ICM are not. The species contrast that anchors the transition
#   claim.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   cross-species modality (single order parameter)
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/theory_crossspecies.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Cross-species figure: the near-critical (unimodal) ICM is conserved in
PRIMATES (human + marmoset) but the mouse ICM is switch-like (bimodal).
Pseudotime-free discriminating test: modality of ICM NANOG-GATA6 (phi).
Human = cESFW atlas 'ICM'; Marmoset = Boroviak ICM_CS3; Mouse = Nowotschin E3.5 inner."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
z = lambda v:(v-v.mean())/(v.std()+1e-9)
def bic(x):
    X=x.reshape(-1,1)
    return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)

# ---- HUMAN ----
LIN=['ICM/TE Branch','ICM','Epi/Hypo Branch','preIm-Epi','Hypo']
SI,C,cG,E,*_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs['a']=ann[m]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); al=pd.DataFrame(ad.X,index=ad.obs_names,columns=ad.var_names)
phiH=(z(al['NANOG'].values)-z(al['GATA6'].values))[ad.obs['a'].values=='ICM']

# ---- MARMOSET (Boroviak ICM_CS3) ----
D='data/nhp/marmoset/RawCounts.csv'; meta=pd.read_csv(D,index_col=0,nrows=2); lineage=meta.iloc[1]
counts=pd.read_csv(D,skiprows=[1,2],index_col=0); counts=counts[~counts.index.duplicated()]
cells=[c for c in counts.columns if lineage.get(c) in ('ICM_CS3','Epi_CS3','Hyp_CS3')]
adm=sc.AnnData(X=counts[cells].T.values.astype(float)); adm.var_names=counts.index; adm.obs_names=cells
adm.obs['lin']=[lineage[c] for c in cells]
sc.pp.normalize_total(adm,target_sum=1e4); sc.pp.log1p(adm); alm=pd.DataFrame(adm.X,index=adm.obs_names,columns=adm.var_names)
phiM=(z(alm['NANOG'].values)-z(alm['GATA6'].values))[np.array(adm.obs['lin'])=='ICM_CS3']

# ---- MOUSE (Nowotschin E3.5 inner) ----
DATA='data/mouse_dev'
F={'E3.5':['GSM3494334_Lib1-3_E3.5_counts.csv.gz','GSM3494335_Lib1-4_E3.5_counts.csv.gz'],
   'E4.5':['GSM3494336_Lib1-1_E4.5_counts.csv.gz','GSM3494337_Lib1-2_E4.5_counts.csv.gz']}
fr,st=[],[]
for s,fs in F.items():
    for f in fs:
        d=pd.read_csv(os.path.join(DATA,f),index_col=0); d.index=[f'{f[:10]}_{i}' for i in range(len(d))]; fr.append(d); st+=[s]*len(d)
mt=pd.concat(fr).fillna(0); adk=sc.AnnData(X=mt.values.astype(float)); adk.var_names=[g.upper() for g in mt.columns]
adk.var_names_make_unique(); adk.obs_names=mt.index; adk.obs['stage']=st
sc.pp.filter_cells(adk,min_genes=1500); sc.pp.normalize_total(adk,target_sum=1e4); sc.pp.log1p(adk)
alk=pd.DataFrame(adk.X,index=adk.obs_names,columns=adk.var_names)
sc.tl.score_genes(adk,[g for g in ['CDX2','GATA3','KRT8','KRT18','TEAD4','EOMES'] if g in adk.var_names],score_name='TE')
inner=(adk.obs['TE'].values<0.15)&(adk.obs['stage'].values=='E3.5')
phiK=(z(alk['NANOG'].values)-z(alk['GATA6'].values))[inner]

data=[("Human ICM","#C0392B",phiH),("Marmoset ICM","#E8873A",phiM),("Mouse ICM (E3.5)","#2E86C1",phiK)]
OUT=os.path.join("plots","revision_results"); os.makedirs(OUT,exist_ok=True)
fig,ax=plt.subplots(1,4,figsize=(15,3.9),gridspec_kw={'width_ratios':[1,1,1,0.9],'wspace':0.32})
for i,(nm,c,x) in enumerate(data):
    b=bic(x); g=np.linspace(x.min()-.7,x.max()+.7,200); d=gaussian_kde(x)(g)
    ax[i].hist(x,bins=18,density=True,color=c,alpha=0.35,edgecolor='white',linewidth=0.5)
    ax[i].plot(g,d,color=c,lw=2.4)
    verdict=("unimodal (single population)" if b<-5 else "bimodal (resolved / switch-like)" if b>5 else "borderline (unresolved; mixture unclear)")
    ax[i].set_title(f"{nm}  (n={len(x)})\n{verdict}  ($\\Delta$BIC={b:+.0f})",fontsize=10)
    ax[i].set_xlabel("$\\phi$ = z(NANOG) - z(GATA6)"); ax[i].set_ylabel("density" if i==0 else "")
    ax[i].axvline(0,color='#999',ls=':',lw=1); ax[i].spines[['top','right']].set_visible(False)
# summary bar
names=[d[0].split(' ICM')[0] for d in data]; bics=[bic(d[2]) for d in data]; cols=[d[1] for d in data]
ax[3].bar(names,bics,color=cols,alpha=0.9,edgecolor='white')
ax[3].axhline(0,color='#333',lw=1); ax[3].set_ylabel("$\\Delta$BIC (>0 bimodal / switch)")
ax[3].set_title("ICM modality across species",fontsize=10)
ax[3].tick_params(axis='x',rotation=20); ax[3].spines[['top','right']].set_visible(False)
ax[3].text(0.5,-0.03,"mouse decisively resolved; primates unresolved (human borderline)",transform=ax[3].transAxes,ha='center',va='top',fontsize=7.5,color='#666')
fig.suptitle("The mouse ICM is resolved into two populations; the primate ICM is a single unresolved population\n"
             "(NANOG-GATA6 order parameter; mouse decisively bimodal, marmoset unimodal, human at the boundary)",
             fontsize=12,y=1.08)
out=os.path.join(OUT,"Fig_crossspecies_criticality.png"); fig.savefig(out,dpi=200,bbox_inches='tight')
print("Saved:",out)
for nm,c,x in data: print(f"  {nm}: n={len(x)} var={x.var():.2f} bic={bic(x):+.2f}")
