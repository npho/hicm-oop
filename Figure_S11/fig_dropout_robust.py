# ============================================================================
# Figure S11  --  dropout bias test + synthetic benchmark
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The single-population reading is not an artefact of the 16% NANOG
#   dropout (dropout is biased, and a synthetic benchmark recovers ground
#   truth).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   dropout bias test + synthetic benchmark
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S11. The human ICM single-population reading is not an artefact
#   of NANOG dropout. (A) Within the ICM, cells with undetected NANOG (n=13)
#   sit at lower ϕ3(multigene composite) than NANOG-detected cells (means
#   -0.47vs +0.01; Mann–Whitney p=0.045), i.e. they are hypoblast-biased:
#   dropout tracks a biological resolution gradient rather than random
#   technical loss. (B) Synthetic benchmark: applying realistic 16% NANOG
#   dropout to simulated ICM populations of the observed size, a genuinely
#   bimodal (two-committed-population) ICM is still recovered as bimodal in
#   100% of simulations and a genuinely unimodal one as unimodal in 99%; the
#   observed human two-gene value (ΔBIC=-2.7) falls in the unimodal regime.
#   Dropout therefore does not mask bimodality, so the single-population
#   classification is not a dropout artefact.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_dropout_robust.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Fig S13: the human ICM single-population reading is not a NANOG-dropout artefact.
(A) NANOG-undetected ICM cells sit at lower phi3 (hypoblast-biased) -> dropout tracks a
    biological gradient, consistent with gradual resolution, not random technical noise.
(B) Synthetic benchmark: with realistic 16% NANOG dropout, a TRULY bimodal ICM is still
    recovered as bimodal; a truly unimodal one as unimodal -> dropout does not hide bimodality."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, sys, os
from sklearn.mixture import GaussianMixture
from scipy import stats
import matplotlib.pyplot as plt
sys.path.insert(0,os.getcwd()); from _load_data import load_embryo_data
z=lambda v:(v-v.mean())/(v.std()+1e-9)
EPI=['NANOG','SOX2','POU5F1','KLF17','TDGF1','DPPA5','PRDM14']; HYP=['GATA6','GATA4','SOX17','PDGFRA','FOXA2','HNF1B']
def dbic(x):
    X=np.asarray(x).reshape(-1,1); return GaussianMixture(1,n_init=3,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=3,random_state=0).fit(X).bic(X)
LIN=['ICM/TE Branch','ICM','Epi/Hypo Branch','preIm-Epi','Hypo']
SI,C,cG,E,*_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values; m=np.isin(ann,LIN)
ad=sc.AnnData(X=C.values[m].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[m]; ad.obs['a']=ann[m]
raw=C.values[m][:, list(C.columns).index('NANOG')]
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
sc.tl.score_genes(ad,[g for g in EPI if g in ad.var_names],score_name='_e'); sc.tl.score_genes(ad,[g for g in HYP if g in ad.var_names],score_name='_h')
phi3=z(ad.obs['_e'].values)-z(ad.obs['_h'].values); icm=ad.obs['a'].values=='ICM'
nz=raw[icm]==0; p3=phi3[icm]
u,pu=stats.mannwhitneyu(p3[nz],p3[~nz],alternative="less")
rng=np.random.RandomState(0); n=int(icm.sum())
def sim(bi):
    if bi:
        k=n//2; na=np.r_[rng.normal(2.5,0.6,k),rng.normal(0.3,0.3,n-k)]; ga=np.r_[rng.normal(0.3,0.3,k),rng.normal(2.5,0.6,n-k)]
    else:
        na=rng.normal(1.2,0.7,n); ga=rng.normal(1.2,0.7,n)
    na=np.clip(na,0,None); ga=np.clip(ga,0,None); d=rng.rand(n)<0.16; na=na.copy(); na[d]=0
    return z(na)-z(ga)
bB=[dbic(sim(1)) for _ in range(300)]; bU=[dbic(sim(0)) for _ in range(300)]

fig,ax=plt.subplots(1,2,figsize=(11,4.3),gridspec_kw={'wspace':0.30})
# A
axA=ax[0]
axA.hist(p3[~nz],bins=14,density=True,alpha=0.55,color='#6C3483',label='NANOG detected (n=%d)'%(~nz).sum())
axA.hist(p3[nz],bins=8,density=True,alpha=0.55,color='#C0392B',label='NANOG undetected (n=%d)'%nz.sum())
axA.axvline(p3[~nz].mean(),color='#6C3483',ls='--'); axA.axvline(p3[nz].mean(),color='#C0392B',ls='--')
axA.set_xlabel('$\\phi_3$ (multigene epiblast $-$ hypoblast)'); axA.set_ylabel('density')
axA.set_title('A  NANOG-undetected ICM cells are hypoblast-biased\n(means %.2f vs %.2f; Mann-Whitney p=%.3f)'%(p3[nz].mean(),p3[~nz].mean(),pu),fontsize=9.8)
axA.legend(fontsize=8,frameon=False); axA.spines[['top','right']].set_visible(False)
# B
axB=ax[1]
axB.hist(bU,bins=30,alpha=0.6,color='#2E86C1',label='true unimodal (transition)')
axB.hist(bB,bins=30,alpha=0.6,color='#E8873A',label='true bimodal (mixture)')
axB.axvline(0,color='#333',ls=':'); axB.axvline(-2.7,color='#6C3483',ls='--',label='observed human ICM ($\\phi_1=-2.7$)')
axB.set_xlabel('$\\Delta$BIC after 16% NANOG dropout ( >0 = bimodal )'); axB.set_ylabel('simulations')
axB.set_title('B  Dropout does not hide bimodality\n(true bimodal recovered as bimodal in %.0f%% of sims)'%(100*np.mean(np.array(bB)>0)),fontsize=9.8)
axB.legend(fontsize=8,frameon=False); axB.spines[['top','right']].set_visible(False)
fig.suptitle('The human ICM single-population reading is not an artefact of NANOG dropout',fontsize=12,y=1.02,fontweight='bold')
out='plots/revision_results/Fig_dropout_robust.png'; fig.savefig(out,dpi=200,bbox_inches='tight',facecolor='white')
print('Saved',out); print('A: undetected mean %.2f vs detected %.2f, p=%.3f'%(p3[nz].mean(),p3[~nz].mean(),pu))
print('B: bimodal recovered %.0f%%, unimodal %.0f%%'%(100*np.mean(np.array(bB)>0),100*np.mean(np.array(bU)<0)))
