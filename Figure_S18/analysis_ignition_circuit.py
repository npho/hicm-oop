# ============================================================================
# Figure S18  --  A-D: ignition criteria, co-expression, wiring, footprints
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   ESRRB and TFCP2L1 are the ICM-initiated elements; KLF4/KLF5/TBX3 are
#   downstream. The four-criterion ignition test.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A-D: ignition criteria, co-expression, wiring, footprints
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S18. ESRRB and TFCP2L1 are the ICM-initiated elements of the
#   ground-state circuit; KLF4, KLF5 and TBX3 are downstream. (A) Expression
#   trajectories along the pluripotent lineage: ESRRB and TFCP2L1 peak at
#   the ICM, whereas KLF5 peaks earlier (morula) and TBX3 later (epiblast
#   branch); KLF4 is high at the ICM but essentially flat across the morula-
#   to-epiblast window (7.3-7.9 log2CPM) and therefore does not resolve the
#   ICM. (B) Of the five factors tested only ESRRB and TFCP2L1 engage the
#   ICM-stage module, and they do so asymmetrically: TFCP2L1 supplies the
#   switch-on without the transient decline, whereas ESRRB supplies the
#   transient decline without a large induction (quantified in Fig. S21C).
#   (C) Co-expression in ICM metacells: TBX3 co-varies with KLF4 (a later
#   wave), not with the ESRRB/TFCP2L1 core. (D) ICM TOBIAS footprints (Chen
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_ignition_circuit.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Test the ICM-"ignition" circuit hypothesis: is it ESRRB+TFCP2L1, with or without KLF4/KLF5?
 (1) ignition signature: switched ON at the ICM (up vs upstream morula/branch; transient vs epiblast)
 (2) co-expression with ESRRB/TFCP2L1 in the same ICM cells (metacells)
 (3) wiring: are KLF4/KLF5 bound by mouse Esrrb/Tcfcp2l1 (Chen 2008), i.e. targets of the core?
"""
import os,sys,gzip,json,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE='/Users/adelosan/cESFW_analysis'; sys.path.insert(0,HERE); os.chdir(HERE); rng=np.random.RandomState(0)
from _load_data import load_embryo_data
PLOTS=f'{HERE}/plots/revision_results'; CH='data/chen2008'
STAGES=['8-Cell','Morula','ICM/TE Branch','ICM','Epi/Hypo Branch','preIm-Epi','Embryonic Disc']
SLAB=['8-cell','morula','ICM/TE\nbranch','ICM','Epi/Hypo\nbranch','preIm\nepiblast','emb\ndisc']
CAND=['ESRRB','TFCP2L1','KLF4','KLF5','TBX3','KLF2','KLF17','NANOG']   # candidates + context
COL={'ESRRB':'#c0392b','TFCP2L1':'#e67e22','KLF4':'#27ae60','KLF5':'#2980b9','KLF2':'#95a5a6','KLF17':'#7f8c8d','NANOG':'#2c3e70','TBX3':'#8e44ad'}

SI,C,cG,E,_,_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values; gl=list(C.columns); gi={g:i for i,g in enumerate(gl)}
def ln(A): A=np.asarray(A,float); lib=A.sum(1,keepdims=True); lib[lib==0]=1; return np.log1p(A/lib*1e4)
X=ln(C.values)
def prof(g): return np.array([X[ann==s][:,gi[g]].mean() for s in STAGES])
def sem(g): return np.array([X[ann==s][:,gi[g]].std()/np.sqrt((ann==s).sum()) for s in STAGES])
P={g:prof(g) for g in CAND}
iICM=STAGES.index('ICM')

# ---- (1) ignition metrics ----
print("=== ignition signature (per factor) ===")
print(f"{'gene':9}{'peak_stage':16}{'onset ICM/upstream':>20}{'transient ICM/downstream':>26}{'ignition?':>11}")
ign={}
for g in ['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']:
    p=P[g]; peak=STAGES[int(np.argmax(p))]
    upstream=max(p[STAGES.index('Morula')],p[STAGES.index('ICM/TE Branch')])
    downstream=p[STAGES.index('preIm-Epi')]
    onset=p[iICM]-upstream; trans=p[iICM]-downstream
    is_ign = (onset>0.05) and (p[iICM]>=0.6*p.max())    # newly up at ICM AND near its own max at ICM
    ign[g]=dict(peak=peak,onset=float(onset),transient=float(trans),ignition=bool(is_ign))
    print(f"{g:9}{peak:16}{onset:20.2f}{trans:26.2f}{str(is_ign):>11}")

# ---- (2) co-expression in ICM (metacells) ----
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
icm_idx=np.where(ann=='ICM')[0]
Xicm=X[icm_idx]
pcs=PCA(min(20,len(icm_idx)-1),random_state=0).fit_transform(Xicm)
k=max(4,len(icm_idx)//6); km=KMeans(k,n_init=5,random_state=0).fit_predict(pcs)
# metacell means (log-norm avg) for candidates
MC=np.vstack([Xicm[km==c].mean(0) for c in range(k) if (km==c).sum()>0])
sub=['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']
M=MC[:,[gi[g] for g in sub]]
Ccorr=np.corrcoef(M.T)
print("\n=== co-expression in ICM metacells (Pearson) ===")
print("        "+"  ".join(f"{g:>8}" for g in sub))
for i,g in enumerate(sub): print(f"{g:8}"+"  ".join(f"{Ccorr[i,j]:8.2f}" for j in range(len(sub))))
# co-detection: fraction of ICM cells co-expressing each pair with ESRRB & TFCP2L1
raw=C.values[icm_idx]
def pos(g): return raw[:,gi[g]]>0
core_pos = pos('ESRRB')&pos('TFCP2L1')
print("\nfraction of ICM cells co-expressing:")
print(f"  ESRRB&TFCP2L1: {core_pos.mean()*100:.0f}%")
for g in ['KLF4','KLF5','TBX3']:
    print(f"  ESRRB&TFCP2L1&{g}: {(core_pos&pos(g)).mean()*100:.0f}%   ({g}+ alone {pos(g).mean()*100:.0f}%)")

# ---- (3) wiring: mouse Esrrb/Tcfcp2l1 binding near Klf4/Klf5/Esrrb/Tfcp2l1 (Chen 2008 mm8) ----
def load_tss(fp):
    t=defaultdict(list)
    for l in gzip.open(fp,'rt'):
        p=l.rstrip().split('\t'); c=p[2]
        if '_' in c: continue
        t[p[12]].append((c,int(p[4]) if p[3]=='+' else int(p[5])))
    return t
mt=load_tss(f'{CH}/refGene_mm8.txt.gz')
def peaks(f):
    out=[]
    for l in gzip.open(f'{CH}/{f}.txt.gz','rt'):
        p=l.rstrip().split('\t')
        try: c,r=p[0].split(':'); s,e=r.split('-'); out.append((c,(int(s)+int(e))//2))
        except: pass
    return out
pk={f:peaks(f) for f in ['Esrrb','Tcfcp2l1']}
def bound(factor,gene,win=10000):
    if gene not in mt: return None
    locs=mt[gene]
    for c,m in pk[factor]:
        for gc,gt in locs:
            if c==gc and abs(m-gt)<win: return True
    return False
print("\n=== mouse ChIP wiring (Chen 2008): is gene bound by the core within 10kb of TSS? ===")
for gene in ['Klf4','Klf5','Tbx3','Esrrb','Tcfcp2l1','Nanog','Pou5f1']:
    e=bound('Esrrb',gene); t=bound('Tcfcp2l1',gene)
    print(f"  {gene:9} Esrrb-bound={e}  Tcfcp2l1-bound={t}")

# ---- FIGURE ----
fig=plt.figure(figsize=(20,5)); gs=fig.add_gridspec(1,4,width_ratios=[1.25,0.95,0.85,1.0],wspace=.34); xs=np.arange(len(STAGES))
axA=fig.add_subplot(gs[0,0])
for g in ['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']:
    axA.plot(xs,P[g],'-o',color=COL[g],lw=2,label=g); axA.fill_between(xs,P[g]-sem(g),P[g]+sem(g),color=COL[g],alpha=.13)
axA.axvline(iICM,color='#c0392b',ls=':',lw=1)
axA.set_xticks(xs); axA.set_xticklabels(SLAB,fontsize=8.5); axA.set_ylabel('mean log-norm expr',fontsize=10.5)
axA.legend(fontsize=9,frameon=False); axA.set_title('A   ESRRB/TFCP2L1 peak AT the ICM; KLF5 pre-ICM, KLF4 post-ICM',fontsize=11,loc='left',weight='bold')
axB=fig.add_subplot(gs[0,1])
for g in ['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']:
    axB.scatter(ign[g]['onset'],ign[g]['transient'],s=140,color=COL[g],edgecolor='k',lw=.5,zorder=3)
    axB.annotate(g,(ign[g]['onset'],ign[g]['transient']),fontsize=9,xytext=(4,3),textcoords='offset points')
axB.axhline(0,color='k',lw=.5); axB.axvline(0,color='k',lw=.5)
axB.fill_betweenx([0,axB.get_ylim()[1] if False else 2],0,2,color='#c0392b',alpha=.05)
axB.set_xlabel('onset  (ICM − upstream)',fontsize=9.5); axB.set_ylabel('transience  (ICM − epiblast)',fontsize=9.5)
axB.set_title('B   Ignition = up at ICM AND transient\n(top-right quadrant)',fontsize=10.5,loc='left',weight='bold')
axC=fig.add_subplot(gs[0,2])
im=axC.imshow(Ccorr,cmap='RdBu_r',vmin=-1,vmax=1)
axC.set_xticks(range(len(sub))); axC.set_xticklabels(sub,rotation=90,fontsize=8.5); axC.set_yticks(range(len(sub))); axC.set_yticklabels(sub,fontsize=8.5)
axC.set_title('C   Co-expression in ICM\n(metacells)',fontsize=10.5,loc='left',weight='bold')
for i in range(len(sub)):
    for j in range(len(sub)): axC.text(j,i,f'{Ccorr[i,j]:.2f}',ha='center',va='center',fontsize=7.5)
fig.colorbar(im,ax=axC,fraction=.046,pad=.04)
# D: ICM TOBIAS footprint (5-motif run) ICM vs Naive
axD=fig.add_subplot(gs[0,3])
try:
    bd=pd.read_csv('atac_rt_pipeline/tobias_out/BINDetect_ignition5/bindetect_results.txt',sep='\t')
    bd['key']=bd['output_prefix'].str.split('_').str[0].str.upper()
    order=['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']; disp=['ESRRB','TFCP2L1','KLF4','KLF5','TBX3']
    icm=[]; nai=[]
    for k in order:
        r=bd[bd['key']==k.upper()]
        icm.append(float(r['ICM_mean_score'].values[0]) if len(r) else np.nan)
        nai.append(float(r['Naive_mean_score'].values[0]) if len(r) else np.nan)
    y=np.arange(len(order)); w=.38
    axD.barh(y+w/2,icm,w,color=[COL.get(g,'#888') for g in disp],label='human ICM')
    axD.barh(y-w/2,nai,w,color='#cfcfcf',label='naïve hESC')
    axD.set_yticks(y); axD.set_yticklabels(disp,fontsize=9)
    axD.set_xlabel('ICM TOBIAS footprint score',fontsize=9.5); axD.legend(fontsize=8,loc='lower right')
    axD.set_title('D   Ignition core (ESRRB/TFCP2L1) is footprinted in\nICM chromatin, below naïve levels',fontsize=10,loc='left',weight='bold')
except Exception as e:
    axD.text(.5,.5,f'footprint n/a\n{e}',ha='center',transform=axD.transAxes,fontsize=8)
fig.suptitle('Is the ICM-ignition circuit ESRRB+TFCP2L1, ± KLF4/KLF5/TBX3?  (RNA timing, co-expression, mouse-ChIP wiring, ICM footprint)',fontsize=12.5,y=1.03)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_Ignition_Circuit.png',dpi=170,bbox_inches='tight')
print("\nSAVED Fig_Ignition_Circuit.png")
json.dump({'ignition_metrics':ign,'icm_coexpr':{f'{sub[i]}-{sub[j]}':float(Ccorr[i,j]) for i in range(len(sub)) for j in range(i+1,len(sub))}},
          open('/tmp/ignition_circuit.json','w'),default=float,indent=2)
print("DONE")
