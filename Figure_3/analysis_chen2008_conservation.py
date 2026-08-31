# ============================================================================
# Figure 3  --  C,D: Chen 2008 mouse ChIP targets mapped to human orthologues
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM ignites a conserved ground-state circuit (ESRRB/TFCP2L1) on a
#   maintained OCT4/SOX2/NANOG core, with mouse ChIP conservation and ICM
#   footprints.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   C,D: Chen 2008 mouse ChIP targets mapped to human orthologues
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_chen2008_conservation.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Conservation of ground-state CIRCUITRY (improved): mouse-ESC ChIP-seq (Chen 2008, GSE11431)
regulons -> human ICM (2PN ATAC + scRNA). Improvements vs first pass:
 - HIGH-CONFIDENCE targets: top-scored, promoter-proximal (<=5kb) peaks only (cut promiscuity)
 - CORRECT contrast: pluripotency (ICM+epiblast) vs TE/differentiated (ICM & epi SHARE the regulon)
 - conserved-core gene list, annotated for ground-state factors
 - integrate w/ RNA three-wave: is the conserved circuit INITIATED at the ICM?
"""
import os,sys,gzip,json,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from collections import defaultdict
from scipy.stats import hypergeom
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
HERE='/Users/adelosan/cESFW_analysis'; sys.path.insert(0,HERE); os.chdir(HERE)
from _load_data import load_embryo_data
CH='data/chen2008'; PLOTS=f'{HERE}/plots/revision_results'; WIN=5000; TOPFRAC=0.5
CORE=['Esrrb','Tcfcp2l1']; CORE_EXT=['Nanog','Oct4','Sox2','Klf4']; MYC=['cMyc','nMyc','E2f1','Zfx']
FACTORS=CORE+CORE_EXT+MYC

def load_tss(fp,c_i,s_i,ts_i,te_i,n_i):
    tss=defaultdict(list)
    with gzip.open(fp,'rt') as f:
        for line in f:
            p=line.rstrip('\n').split('\t'); c=p[c_i]
            if '_' in c: continue
            t=int(p[ts_i]) if p[s_i]=='+' else int(p[te_i]); tss[p[n_i]].append((c,t))
    return tss
def chrom_index(tss):
    ct=defaultdict(list)
    for g,locs in tss.items():
        for c,t in locs: ct[c].append((t,g))
    pos={}; sym={}
    for c in ct:
        v=sorted(ct[c]); pos[c]=np.array([x[0] for x in v]); sym[c]=[x[1] for x in v]
    return pos,sym
mtss=load_tss(f'{CH}/refGene_mm8.txt.gz',2,3,4,5,12); mpos,msym=chrom_index(mtss)

def targets(peakfile,win=WIN,top_frac=TOPFRAC):
    peaks=[]
    with gzip.open(peakfile,'rt') as f:
        for line in f:
            p=line.rstrip('\n').split('\t'); loc=p[0]
            try:
                c,rng=loc.split(':'); s,e=rng.split('-'); mid=(int(s)+int(e))//2; sc=float(p[3]) if len(p)>3 else 1.0
            except: continue
            if c in mpos: peaks.append((c,mid,sc))
    if not peaks: return set()
    peaks.sort(key=lambda x:-x[2]); keep=peaks[:max(1,int(len(peaks)*top_frac))]   # top-scored
    tg=set()
    for c,mid,_ in keep:
        lo=np.searchsorted(mpos[c],mid-win); hi=np.searchsorted(mpos[c],mid+win)
        for k in range(lo,hi): tg.add(msym[c][k])
    return tg
mouse_tg={f:targets(f'{CH}/{f}.txt.gz') for f in FACTORS}
for f in FACTORS: print(f"{f:9} hi-conf proximal targets {len(mouse_tg[f])}")

SI,C,cG,E,_,_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values
hgenes=set(C.columns); glist=list(C.columns)
def ln(A): A=np.asarray(A,float); lib=A.sum(1,keepdims=True); lib[lib==0]=1; return np.log1p(A/lib*1e4)
X=ln(C.values)
def to_human(mg): return {g.upper() for g in mg}&hgenes
hum_tg={f:to_human(mouse_tg[f]) for f in FACTORS}

# programs: pluripotency (ICM+epiblast) vs TE (correct contrast)
icm=X[ann=='ICM'].mean(0); epi=X[ann=='preIm-Epi'].mean(0); te=X[np.isin(ann,['Mid TE','Polar TE','Mural TE'])].mean(0)
pluri=(icm+epi)/2
def upset(a,b,n=400):
    d=a-b; ex=np.maximum(a,b)>0.3; g=np.array(glist)[ex]; return set(g[np.argsort(-d[ex])][:n])
PLURI_up=upset(pluri,te); DIFF_up=upset(te,pluri)
ICM_active=set(np.array(glist)[icm>0.5])
print(f"PLURI-up {len(PLURI_up)} | DIFF(TE)-up {len(DIFF_up)}")

# human 2PN ICM ATAC accessible genes
htss=load_tss('data/GSE101571/refGene_hg19.txt.gz',2,3,4,5,12); hpos,hsym=chrom_index(htss)
acc=set()
with gzip.open('data/gao2018_atac/icm_2pn_peaks.bed.gz','rt') as f:
    for line in f:
        p=line.split('\t'); c=p[0]
        if c not in hpos: continue
        mid=(int(p[1])+int(p[2]))//2; lo=np.searchsorted(hpos[c],mid-WIN); hi=np.searchsorted(hpos[c],mid+WIN)
        for k in range(lo,hi): acc.add(hsym[c][k])
acc&=hgenes

UNI=hgenes; N=len(UNI)
def enr(t,s):
    t=t&UNI; s=s&UNI; k=len(t&s)
    if not t or not s: return 0.0,k,len(t)
    exp=len(t)*len(s)/N; return (k+.5)/(exp+.5),k,len(t)
rows=[]
for f in FACTORS:
    o_pl=enr(hum_tg[f],PLURI_up)[0]; o_df=enr(hum_tg[f],DIFF_up)[0]; o_ac=enr(hum_tg[f],acc)[0]
    mod='Core' if f in CORE else ('Core-ext' if f in CORE_EXT else 'Myc')
    rows.append(dict(factor=f,module=mod,n=enr(hum_tg[f],acc)[2],pluri_or=o_pl,diff_or=o_df,acc_or=o_ac,pluri_directed=o_pl/max(o_df,1e-9)))
R=pd.DataFrame(rows); print(R[['factor','module','n','pluri_or','diff_or','pluri_directed','acc_or']].round(2).to_string())

# conserved core: Esrrb & Tcfcp2l1 co-targets, pluri-program + ICM-accessible
core_co=hum_tg['Esrrb']&hum_tg['Tcfcp2l1']
conserved=sorted(core_co & (ICM_active|PLURI_up) & acc)
GSGENES=['ESRRB','TFCP2L1','KLF2','KLF4','KLF5','TBX3','PRDM14','ZFP42','NR0B1','DPPA3','DPPA5','KLF17','SALL4','NANOG','POU5F1','SOX2','GBX2','FGF4','TDGF1','GDF3','TFCP2L1']
gs_in=[g for g in conserved if g in GSGENES]
print(f"\nCONSERVED CORE (Esrrb&Tcfcp2l1 co-targets, pluri-active + ICM-accessible in human): {len(conserved)}")
print("ground-state members captured:",gs_in)

# integrate with RNA three-wave: is the conserved circuit INITIATED at ICM?
init_at_icm=[]
stages=['Morula','ICM/TE Branch','ICM','preIm-Epi','Embryonic Disc']
for g in gs_in:
    gi=glist.index(g); prof=[X[ann==s][:,gi].mean() for s in stages]
    init_at_icm.append((g,stages[int(np.argmax(prof))]))
print("peak stage of captured ground-state conserved genes:",init_at_icm)

# ---- FIGURE ----
fig,axs=plt.subplots(1,2,figsize=(15,5.8),gridspec_kw={'width_ratios':[1.1,1]})
CM={'Core':'#c0392b','Core-ext':'#e67e22','Myc':'#2c3e70'}
ax=axs[0]; Rs=R.sort_values('pluri_directed'); y=np.arange(len(Rs))
ax.barh(y,np.log2(Rs['pluri_directed']),color=[CM[m] for m in Rs['module']],edgecolor='k',lw=.4)
ax.set_yticks(y); ax.set_yticklabels(Rs['factor'],fontsize=9.5); ax.axvline(0,color='k',lw=.7)
ax.set_xlabel('pluripotency-directedness of targets  log2( PLURI-up OR / TE-up OR )',fontsize=9.5)
ax.set_title('C   Core ground-state regulon (Esrrb/Tcfcp2l1 + Oct4/Sox2/Nanog) targets are\nPLURIPOTENCY-directed; Myc targets are not',fontsize=10.5,loc='left',weight='bold')
ax.legend(handles=[Patch(color=c,label=m) for m,c in CM.items()],fontsize=9,loc='lower right')
ax=axs[1]
for _,r in R.iterrows():
    ax.scatter(np.log2(r['pluri_or']),np.log2(r['acc_or']),s=130,color=CM[r['module']],edgecolor='k',lw=.5,zorder=3)
    ax.annotate(r['factor'],(np.log2(r['pluri_or']),np.log2(r['acc_or'])),fontsize=8.5,xytext=(4,3),textcoords='offset points')
ax.axhline(0,color='k',lw=.5); ax.axvline(0,color='k',lw=.5)
ax.set_xlabel('target enrichment in human PLURI program  log2(OR)',fontsize=9.5)
ax.set_ylabel('target enrichment in ICM-ACCESSIBLE (2PN) genes  log2(OR)',fontsize=9.5)
ax.set_title(f'D   Two-layer conservation. {len(conserved)} conserved co-target genes;\nground-state members: {", ".join(gs_in[:8])}',fontsize=10,loc='left',weight='bold')
fig.suptitle('Conserved ground-state circuitry: mouse-ESC ChIP (Chen 2008, GSE11431, hi-conf) -> human pluripotent compartment / ICM',fontsize=11,y=1.02)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_Chen2008_Conservation.png',dpi=175,bbox_inches='tight')
print("SAVED Fig_Chen2008_Conservation.png")
json.dump({'window':WIN,'top_frac':TOPFRAC,'enrichment':R.to_dict('records'),
           'n_conserved':len(conserved),'ground_state_conserved':gs_in,'init_stage':init_at_icm,
           'conserved_examples':conserved[:80]},open('/tmp/chen2008_conservation.json','w'),default=float,indent=2)
print("DONE")
