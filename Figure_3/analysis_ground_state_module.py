# ============================================================================
# Figure 3  --  A,B: ground-state module across the lineage
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM ignites a conserved ground-state circuit (ESRRB/TFCP2L1) on a
#   maintained OCT4/SOX2/NANOG core, with mouse ChIP conservation and ICM
#   footprints.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A,B: ground-state module across the lineage
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 3. The human ICM initiates a conserved pluripotency circuit.
#   (A,B) A conserved murine ground-state module in the human ICM: the
#   ESRRB/TFCP2L1 module peaks across the ICM on a maintained
#   OCT4/SOX2/NANOG core; TFCP2L1 supplies the switch-on (+2.8 log2CPM from
#   the morula) and is then maintained, whereas ESRRB is already present in
#   the morula, rises only modestly (+0.6) and is instead distinguished by
#   its extinction at the epiblast (reproduced across three studies, Fig.
#   S19; quantified in Fig. S21C). (C,D) Mouse ESC ESRRB/TFCP2L1 targets
#   (Chen et al., 2008) mapped to human orthologues define a conserved
#   pluripotency regulon distinct from a MYC control module; about half of
#   the co-targets are conserved. The 354 conserved co-targets are dominated
#   by broadly expressed genes but include recognizable pluripotency and
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_ground_state_module.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Ground-state regulatory MODULE analysis (RNA-only), on the human embryo atlas.
 GS1  (A) CORE vs GROUND-STATE arm activity across the pluripotency lineage (AUCell, metacells)
      (B) per-gene z-scored remodeling heatmap
      (C) canonical ground-state TF trajectories (ESRRB, TFCP2L1) vs maintained core (NANOG)
 GS2  module COHERENCE per stage vs detection-matched null ("preservation test"),
      single-cell (power) + TE negative contrast, + human ICM co-expression heatmap
HONEST framing: the CORE pluripotency network (OCT4/SOX2/NANOG/SALL4) is established early and
maintained; the mouse-ground-state-SPECIFIC factors (ESRRB, TFCP2L1) are TRANSIENTLY maximal at
the ICM and remodel toward the epiblast -> "transient deployment of a ground-state architecture".
"""
import os,sys,warnings,json; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE='/Users/adelosan/cESFW_analysis'; sys.path.insert(0,HERE); os.chdir(HERE)
from _load_data import load_embryo_data
PLOTS=f'{HERE}/plots/revision_results'; os.makedirs(PLOTS,exist_ok=True); rng=np.random.RandomState(0)

CORE  =['POU5F1','SOX2','NANOG','SALL4']
GROUND=['ESRRB','TFCP2L1','KLF4','KLF5','TBX3','ZFP42','KLF17','DPPA3','DPPA5']
MODULE=CORE+GROUND
STAGES=['8-Cell','Morula','ICM/TE Branch','ICM','Epi/Hypo Branch','preIm-Epi','Embryonic Disc']
SLAB ={'8-Cell':'8-cell','Morula':'morula','ICM/TE Branch':'ICM/TE\nbranch','ICM':'ICM',
       'Epi/Hypo Branch':'Epi/Hypo\nbranch','preIm-Epi':'preIm\nepiblast','Embryonic Disc':'embryonic\ndisc','Mid TE':'TE\n(contrast)'}
COH_STAGES=STAGES+['Mid TE']

def lognorm(A):
    A=np.asarray(A,float); lib=A.sum(1,keepdims=True); lib[lib==0]=1; return np.log1p(A/lib*1e4)
def aucell(counts,gidx,top_frac=0.05):
    X=np.asarray(counts,float); n=X.shape[1]; thr=int(np.ceil(top_frac*n))
    order=np.argsort(-X,axis=1); ranks=np.empty_like(order); rows=np.arange(X.shape[0])[:,None]
    ranks[rows,order]=np.arange(1,n+1)[None,:]; mr=ranks[:,gidx]; nmod=mr.shape[1]; auc=np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        r=mr[i]; r=r[r<=thr]; auc[i]=np.sum(thr-r+1)
    return auc/(thr*nmod)
def metacells(Xln,labels,target=8,seed=0):
    from sklearn.cluster import KMeans; from sklearn.decomposition import PCA
    idx=[]; lab=[]
    for L in pd.unique(labels):
        ii=np.where(labels==L)[0]; k=max(2,int(round(len(ii)/target)))
        if len(ii)<=k:
            for j in ii: idx.append([j]); lab.append(L)
            continue
        P=PCA(min(20,len(ii)-1),random_state=seed).fit_transform(Xln[ii])
        km=KMeans(k,n_init=5,random_state=seed).fit_predict(P)
        for c in range(k):
            m=ii[km==c]
            if len(m): idx.append(list(m)); lab.append(L)
    return idx,np.array(lab)
def offdiag(M):
    if M.shape[0]<3: return np.nan
    Cc=np.corrcoef(M.T); iu=np.triu_indices(Cc.shape[0],1); return np.nanmean(Cc[iu])

# ---------------- load human ----------------
SI,C,cG,E,_,_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values
genes=list(C.columns); gi={g:k for k,g in enumerate(genes)}
mod=[g for g in MODULE if g in gi]; core=[g for g in CORE if g in gi]; ground=[g for g in GROUND if g in gi]
sub=np.isin(ann,COH_STAGES); Craw=C.values[sub]; a=ann[sub]; Xln=lognorm(Craw)
gidx=[gi[g] for g in mod]; cidx=[gi[g] for g in core]; grdx=[gi[g] for g in ground]
mc_idx,mc_lab=metacells(Xln,a,target=8); MCraw=np.vstack([Craw[m].sum(0) for m in mc_idx]); MCln=lognorm(MCraw)
print("metacells",len(mc_idx))
auc_core=aucell(MCraw,cidx); auc_grnd=aucell(MCraw,grdx)

# per-gene z across STAGES (single-cell means)
gene_stage=np.vstack([Xln[a==s][:,gidx].mean(0) for s in STAGES]); z=(gene_stage-gene_stage.mean(0))/(gene_stage.std(0)+1e-9)

# ------- FIGURE GS1 (A heatmap | B ESRRB/TFCP2L1 trajectory) -------
fig=plt.figure(figsize=(14.5,5.6)); gs=fig.add_gridspec(1,2,width_ratios=[1.25,1.0],wspace=.28); xs=np.arange(len(STAGES))
# A heatmap
axB=fig.add_subplot(gs[0,0]); og=np.argsort(-np.argmax(z,axis=0))
im=axB.imshow(z[:,og].T,aspect='auto',cmap='RdBu_r',vmin=-1.8,vmax=1.8)
axB.set_yticks(range(len(mod))); axB.set_yticklabels([np.array(mod)[og][k] for k in range(len(mod))],fontsize=9)
axB.set_xticks(range(len(STAGES))); axB.set_xticklabels([SLAB[s] for s in STAGES],fontsize=9)
axB.axvline(STAGES.index('ICM'),color='k',ls=':',lw=1)
axB.set_title('A   Per-factor remodeling across the lineage (row z-score)',fontsize=12,loc='left',weight='bold')
fig.colorbar(im,ax=axB,fraction=.03,pad=.02).set_label('z (mean expr)',fontsize=9)
# B trajectories ESRRB TFCP2L1 vs NANOG
axC=fig.add_subplot(gs[0,1])
for g,col,ls in [('ESRRB','#c0392b','-'),('TFCP2L1','#e67e22','-'),('NANOG','#2c3e70','--')]:
    m=np.array([Xln[a==s][:,gi[g]].mean() for s in STAGES]); se=np.array([Xln[a==s][:,gi[g]].std()/np.sqrt((a==s).sum()) for s in STAGES])
    axC.plot(xs,m,ls,color=col,lw=2.3,label=g); axC.fill_between(xs,m-se,m+se,color=col,alpha=.15)
axC.axvline(STAGES.index('ICM'),color='#c0392b',ls=':',lw=1,alpha=.6)
axC.set_xticks(xs); axC.set_xticklabels([SLAB[s] for s in STAGES],fontsize=9)
axC.set_ylabel('mean log-norm expression',fontsize=11); axC.legend(fontsize=10,loc='upper right',frameon=False)
axC.set_title('B   Ground-state TFs ESRRB / TFCP2L1 peak at the ICM\nand decline; NANOG rises to the disc',fontsize=12,loc='left',weight='bold')
fig.suptitle('Curated murine ground-state module in the human embryo (RNA-only): ground-state factors are transiently '
             'deployed at the ICM\non a retained OCT4/SOX2/NANOG core',fontsize=12.5,y=1.06)
fig.savefig(f'{PLOTS}/Fig_GroundState_Module_Lineage.png',dpi=175,bbox_inches='tight'); print("SAVED GS1")

# ---------------- GS2 coherence (single-cell power + TE contrast + null) ----------------
def stage_coherence_sc(stage,nperm=1000):
    ii=np.where(a==stage)[0]; M=Xln[ii][:,gidx]; obs=offdiag(M)
    det=(Craw[ii]>0).mean(0); dec=np.clip((det*10).astype(int),0,9); pools={d:np.where(dec==d)[0] for d in range(10)}
    for d in range(10):
        if len(pools[d])==0: pools[d]=np.arange(Craw.shape[1])
    null=np.array([offdiag(Xln[ii][:,[rng.choice(pools[dec[gi[g]]]) for g in mod]]) for _ in range(nperm)])
    zc=(obs-np.nanmean(null))/(np.nanstd(null)+1e-9)
    return obs,np.nanmean(null),np.nanpercentile(null,2.5),np.nanpercentile(null,97.5),zc,len(ii)
coh={s:stage_coherence_sc(s) for s in COH_STAGES}
for s in COH_STAGES: print(f"coh {s:16} obs={coh[s][0]:.3f} null={coh[s][1]:.3f} z={coh[s][4]:.1f} n={coh[s][5]}")
# NOTE: co-expression coherence is CONFOUNDED (Mid TE shows the highest coherence via an on/off
# gradient, even after depth control) -> NOT presented as a figure; reported honestly in text.
M_icm=MCln[mc_lab=='ICM'][:,gidx]; Ch=np.corrcoef(M_icm.T)
np.save('/tmp/human_icm_modcorr.npy',Ch); np.save('/tmp/human_mod.npy',np.array(mod))
# store per-stage magnitude of ESRRB/TFCP2L1 for cross-species panel
mag={g:{s:float(Xln[a==s][:,gi[g]].mean()) for s in STAGES} for g in mod}
json.dump(mag,open('/tmp/human_mod_magnitude.json','w'))
json.dump({'coherence':{s:{'obs':float(coh[s][0]),'null':float(coh[s][1]),'z':float(coh[s][4]),'n':int(coh[s][5])} for s in COH_STAGES},
           'module':mod,'core':core,'ground':ground,
           'arm_activity':{s:{'core':float(auc_core[mc_lab==s].mean()),'ground':float(auc_grnd[mc_lab==s].mean())} for s in STAGES},
           'ICM_meanr':float(offdiag(M_icm))},open('/tmp/ground_state_human.json','w'),indent=2)
print("DONE human")
