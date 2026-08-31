# ============================================================================
# Figure 2  --  progressive assembly of the naive programme
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The canonical naive network is assembled only AFTER the ICM: ICM is
#   demethylation-permissive (high TET2, near-absent UHRF1) and low MYC.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   progressive assembly of the naive programme
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 2. The canonical human naïve pluripotency network is assembled
#   only after the ICM. Bottom panels: the ICM is demethylation-permissive,
#   with high TET2 and near-absent UHRF1, whereas DNMT3B rises and TET2
#   declines during epiblast formation (robust in Petropoulos and Yanagida,
#   the two studies adequately sampling the resolved epiblast, with
#   directional support from Meistermann; Fig. S19). c-MYC remains low
#   across the preimplantation pluripotent compartment (ICM and epiblast
#   alike) and increases only after implantation. Thus, ground-state
#   epigenetic hallmarks precedes assembly of the canonical naïve
#   transcription-factor network in humans (unlike the mouse, where they
#   coincide). These trajectories are reproducible across three embryo
#   datasets (Fig. S19).
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_progressive_assembly.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Stress-test the "progressive assembly of pluripotency" model against the human atlas:
 (1) two-wave timing: per-gene peak stage for naive-TF / naive-signalling / hypoblast / TE modules
 (2) module scores across morula->ICM->epiblast (predictions A/B/C)
 (3) coherence reality-check: does within-module co-expression INCREASE from ICM to epiblast?
"""
import os,sys,warnings,json; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
HERE='/Users/adelosan/cESFW_analysis'; sys.path.insert(0,HERE); os.chdir(HERE)
from _load_data import load_embryo_data
PLOTS=f'{HERE}/plots/revision_results'; rng=np.random.RandomState(0)

MODULES={
 'core naive TF':['ESRRB','TFCP2L1','KLF2','KLF4','KLF5','NANOG','TBX3','GBX2','SALL4','POU5F1','SOX2'],
 'naive signalling/output':['NODAL','TDGF1','LEFTY1','LEFTY2','DPPA3','DPPA5','DNMT3L','ZFP42','GDF3','KLF17'],
 'hypoblast competence':['GATA6','GATA4','SOX17','PDGFRA'],
 'TE exclusion':['GATA3','TFAP2C','CDX2'],
}
MCOL={'core naive TF':'#c0392b','naive signalling/output':'#e67e22','hypoblast competence':'#2ecc71','TE exclusion':'#7f8c8d'}
STAGES=['8-Cell','Morula','ICM/TE Branch','ICM','Epi/Hypo Branch','preIm-Epi','Embryonic Disc']
SLAB=['8-cell','morula','ICM/TE\nbranch','ICM','Epi/Hypo\nbranch','preIm\nepiblast','embryonic\ndisc']

SI,C,cG,E,_,_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values
genes=list(C.columns); gi={g:k for k,g in enumerate(genes)}
def ln(A): A=np.asarray(A,float); lib=A.sum(1,keepdims=True); lib[lib==0]=1; return np.log1p(A/lib*1e4)
X=ln(C.values)
present={m:[g for g in gs if g in gi] for m,gs in MODULES.items()}
allg=[g for m in present for g in present[m]]; missing=[g for m in MODULES for g in MODULES[m] if g not in gi]
print("missing:",missing)

# per-gene stage means (z across stages)
gm=np.vstack([np.array([X[ann==s][:,gi[g]].mean() for s in STAGES]) for g in allg])
z=(gm-gm.mean(1,keepdims=True))/(gm.std(1,keepdims=True)+1e-9)
gene_mod=[m for m in present for g in present[m]]
peak_stage=np.argmax(z,axis=1)
order=np.argsort(peak_stage+0.01*np.arange(len(allg)))   # by peak stage

# module scores per stage (mean of z)
modscore={m:np.array([z[[i for i,gg in enumerate(gene_mod) if gg==m]][:,s].mean() for s in range(len(STAGES))]) for m in present}

# coherence ICM vs preIm-Epi (naive TF+signalling), single-cell, detection-matched null
naiveset=[g for g in present['core naive TF']+present['naive signalling/output']]; nidx=[gi[g] for g in naiveset]
def offdiag(M):
    Cc=np.corrcoef(M.T); iu=np.triu_indices(Cc.shape[0],1); return np.nanmean(Cc[iu])
def coh(stage,nperm=500):
    ii=np.where(ann==stage)[0]; obs=offdiag(X[ii][:,nidx])
    det=(C.values[ii]>0).mean(0); dec=np.clip((det*10).astype(int),0,9); pools={d:np.where(dec==d)[0] for d in range(10)}
    for d in range(10):
        if len(pools[d])==0: pools[d]=np.arange(len(genes))
    null=np.array([offdiag(X[ii][:,[rng.choice(pools[dec[gi[g]]]) for g in naiveset]]) for _ in range(nperm)])
    return obs,(obs-np.nanmean(null))/(np.nanstd(null)+1e-9),len(ii)
coh_icm=coh('ICM'); coh_epi=coh('preIm-Epi'); coh_mor=coh('Morula')
print(f"coherence Morula obs={coh_mor[0]:.3f} z={coh_mor[1]:.1f} | ICM obs={coh_icm[0]:.3f} z={coh_icm[1]:.1f} | preIm-Epi obs={coh_epi[0]:.3f} z={coh_epi[1]:.1f}")

# ---- FIGURE ----
fig=plt.figure(figsize=(17,6.6)); g=fig.add_gridspec(1,3,width_ratios=[1.35,1.0,0.7],wspace=.3)
# A heatmap ordered by peak stage, module color strip
axA=fig.add_subplot(g[0,0])
im=axA.imshow(z[order],aspect='auto',cmap='RdBu_r',vmin=-1.6,vmax=1.6)
axA.set_yticks(range(len(allg))); axA.set_yticklabels([allg[i] for i in order],fontsize=7.5)
axA.set_xticks(range(len(STAGES))); axA.set_xticklabels(SLAB,fontsize=8.5)
axA.axvline(STAGES.index('ICM'),color='k',ls=':',lw=1)
for k,i in enumerate(order):
    axA.add_patch(plt.Rectangle((-1.0,k-0.5),0.8,1,color=MCOL[gene_mod[i]],clip_on=False))
axA.set_xlim(-1.0,len(STAGES)-0.5)
axA.set_title('A   Genes ordered by peak stage — waves of assembly',fontsize=11.5,loc='left',weight='bold')
axA.legend(handles=[Patch(color=c,label=m) for m,c in MCOL.items()],fontsize=7.5,loc='lower right',framealpha=.9)
fig.colorbar(im,ax=axA,fraction=.03,pad=.02).set_label('z',fontsize=8)
# B module scores
axB=fig.add_subplot(g[0,1]); xs=np.arange(len(STAGES))
for m in present:
    axB.plot(xs,modscore[m],'-o',color=MCOL[m],lw=2.2,label=m)
axB.axvline(STAGES.index('ICM'),color='k',ls=':',lw=1)
axB.set_xticks(xs); axB.set_xticklabels(SLAB,fontsize=8.5); axB.set_ylabel('module score (mean z)',fontsize=10.5)
axB.legend(fontsize=8,loc='lower left'); axB.set_title('B   Signalling completes at epiblast, hypoblast falls;\nbut the whole "naive-TF module" is FLAT (only ESRRB/TFCP2L1 initiate at ICM)',fontsize=10,loc='left',weight='bold')
# C coherence check
axC=fig.add_subplot(g[0,2])
vals=[coh_mor[0],coh_icm[0],coh_epi[0]]; zs=[coh_mor[1],coh_icm[1],coh_epi[1]]
axC.bar(range(3),vals,color=['#8fbfe0','#c0392b','#2c3e70'])
for i,(v,zz) in enumerate(zip(vals,zs)): axC.text(i,v+.001,f'z={zz:.1f}',ha='center',fontsize=9)
axC.set_xticks(range(3)); axC.set_xticklabels(['morula','ICM','preIm\nepiblast'],fontsize=9)
axC.set_ylabel('naive-module coherence',fontsize=10)
axC.set_title('C   Coherence rises modestly morula→ICM→epiblast\n(weak; +confounded — see TE artifact)',fontsize=10,loc='left',weight='bold')
fig.suptitle('Testing the "progressive assembly of pluripotency" model on the human atlas (RNA)',fontsize=13,y=1.02)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_ProgressiveAssembly.png',dpi=170,bbox_inches='tight')
print("SAVED Fig_ProgressiveAssembly.png")

# report per-gene peak stage by module
print("\n== peak stage by gene ==")
for m in present:
    for g_ in present[m]:
        i=allg.index(g_); print(f"  {m:24} {g_:8} peaks at {STAGES[peak_stage[i]]}")
json.dump({'peak_stage':{allg[i]:STAGES[peak_stage[i]] for i in range(len(allg))},
           'coherence':{'morula':coh_mor[:2],'ICM':coh_icm[:2],'preIm-Epi':coh_epi[:2]}},
          open('/tmp/progressive_assembly.json','w'),default=float,indent=2)
