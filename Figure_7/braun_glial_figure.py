# ============================================================================
# Figure 7  --  A,B: glial astrocyte/oligodendrocyte decision
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The transition-versus-state framework generalises to four further human
#   fate decisions, each with a synthetic 50/50 mixture as positive control.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A,B: glial astrocyte/oligodendrocyte decision
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 7. The transition-versus-state framework generalises across four
#   human fate transitions, from development to regeneration. Order-
#   parameter modality test applied to (A,B) the astrocyte/oligodendrocyte
#   choice of glial progenitors in first-trimester human brain (Braun et
#   al., 2023; thalamus), (C,D) the absorptive/secretory choice of
#   intestinal stem cells in adult small intestine (Burclaff et al., 2022),
#   (E,F) the syncytiotrophoblast/extravillous choice of villous
#   cytotrophoblast in the placenta (Arutyunyan et al., 2023), and (G,H) the
#   AT2-to-AT1 conversion of the KRT8-positive alveolar intermediate in
#   fibrotic adult lung (Habermann et al., 2020). (A,C,E,G) φ = z(fate-A) −
#   z(fate-B): in each system the progenitor or intermediate is a single,
#   unimodal population positioned between its two committed fates.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/figure_reconstruction/braun_glial_figure.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import pandas as pd, numpy as np, matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
z=lambda v:(v-v.mean())/(v.std()+1e-9); rng=np.random.RandomState(0)
def d1(x):
    X=np.asarray(x).reshape(-1,1); return GaussianMixture(1,n_init=2,random_state=0).fit(X).bic(X)-GaussianMixture(2,n_init=2,random_state=0).fit(X).bic(X)
def boot(x,n=150,reps=80):
    x=np.asarray(x)
    if len(x)<n: return np.nan
    return np.median([d1(x[rng.choice(len(x),n,replace=False)]) for _ in range(reps)])
E=pd.read_csv('braun_markers.tsv',sep='\t',index_col=0); E.index=[g.split('|')[0] for g in E.index]; E=E.T
m=pd.read_csv('braun_meta.tsv',sep='\t',low_memory=False).set_index('cell_name')
common=E.index.intersection(m.index); E=E.loc[common]; m=m.loc[common]
numi=pd.to_numeric(m['numi'],errors='coerce').fillna(pd.to_numeric(m.get('ncount_rna'),errors='coerce')).astype(float).values; numi[numi<=0]=np.nan
En=pd.DataFrame(np.log1p(E.values/numi[:,None]*1e4),index=E.index,columns=E.columns)
AC=[g for g in ['AQP4','GFAP','AGT','ALDH1L1','GJA1','S100B'] if g in En.columns]; OL=[g for g in ['OLIG1','OLIG2','SOX10','PDGFRA','CSPG4','PLP1','MBP'] if g in En.columns]
m=m.assign(astro=En[AC].mean(1).values, oligo=En[OL].mean(1).values)
g=m[(m['cell_class']=='glia')&(m['structure']=='thalamus')].copy()
phi=z(g['astro'].values)-z(g['oligo'].values); ct=g['cell_type'].values
pops=[('radial glia','glial progenitor','#7D3C98'),('astrocyte','astrocyte (fate A)','#C0392B'),('oligodendrocyte','oligodendrocyte (fate B)','#2E86C1')]
fig,ax=plt.subplots(1,2,figsize=(11,4.3),gridspec_kw={'wspace':0.28,'width_ratios':[1.3,1]})
gr=np.linspace(phi.min()-.5,phi.max()+.5,200); dbics={}
for nm,lab,c in pops:
    x=phi[ct==nm]; dbics[nm]=boot(x)
    ax[0].plot(gr,gaussian_kde(x)(gr),color=c,lw=2.4,label=f'{lab} (n={len(x)})'); ax[0].axvline(np.mean(x),color=c,ls=':',lw=1.2)
ax[0].set_xlabel('$\\phi$ = z(astrocyte) $-$ z(oligodendrocyte)'); ax[0].set_ylabel('density')
ax[0].set_title('Human glial progenitor is a single unresolved\npopulation between the two fates (thalamus)',fontsize=10)
ax[0].legend(fontsize=8,frameon=False); ax[0].spines[['top','right']].set_visible(False)
names=[p[1].split(' (')[0] for p in pops]; vals=[dbics[p[0]] for p in pops]; cols=[p[2] for p in pops]
ax[1].bar(range(3),vals,color=cols,alpha=0.9,edgecolor='white'); ax[1].axhline(0,color='#333',lw=1); ax[1].axhspan(-5,5,color='#ccc',alpha=.35)
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(['glial\nprogenitor','astrocyte','oligo'],fontsize=9)
ax[1].set_ylabel('calibrated $\\Delta$BIC (>0 bimodal)'); ax[1].set_title('Progenitor unimodal (not a resolved mixture);\nn=150 subsample, 80 reps',fontsize=10)
ax[1].spines[['top','right']].set_visible(False)
fig.suptitle('The transition-vs-state test generalises to a human neural bipotent decision (Braun et al. 2023)',fontsize=12,y=1.02,fontweight='bold')
out='/Users/adelosan/cESFW_analysis/plots/revision_results/Fig_braun_glial.png'; fig.savefig(out,dpi=200,bbox_inches='tight',facecolor='white')
print('Saved',out); print('dBIC:',{k:round(v,1) for k,v in dbics.items()}); print('mean phi:',{p[0]:round(float(np.mean(phi[ct==p[0]])),2) for p in pops})
