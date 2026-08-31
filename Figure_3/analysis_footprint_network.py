# ============================================================================
# Figure 3  --  J,K: TOBIAS footprints (NANOG/GATA6/ESRRB)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM ignites a conserved ground-state circuit (ESRRB/TFCP2L1) on a
#   maintained OCT4/SOX2/NANOG core, with mouse ChIP conservation and ICM
#   footprints.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   J,K: TOBIAS footprints (NANOG/GATA6/ESRRB)
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_footprint_network.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Footprint-anchored regulatory analysis (the stronger, defensible 'network' the user asked for):
 1. For each TF, take its TOBIAS ICM-bound footprints -> assign to nearest gene TSS (refGene hg19).
 2. From the RNA atlas, define ICM-up and EPI-up gene programs (ICM vs preimplantation epiblast).
 3. Test which TFs' footprint-target genes are enriched for the ICM program vs the EPI program
    (hypergeometric) -> which regulators drive each side of the ICM->epiblast transition.
 4. Concordance: do sites GAINED in ICM chromatin (ICM>naive footprint) sit near ICM-up genes?
"""
import os,sys,gzip,json,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import hypergeom
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE='/Users/adelosan/cESFW_analysis'; sys.path.insert(0,HERE); os.chdir(HERE)
from _load_data import load_embryo_data
BD='atac_rt_pipeline/tobias_out/BINDetect'; PLOTS=f'{HERE}/plots/revision_results'
WIN=50000   # site-to-TSS window

# ---------- refGene hg19 TSS ----------
tss={}   # symbol -> list of (chrom, tss)
with gzip.open('data/GSE101571/refGene_hg19.txt.gz','rt') as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        # UCSC refGene: 2 name,3 chrom,4 strand,5 txStart,6 txEnd, 13 name2
        chrom=p[2]; strand=p[3]; t=int(p[4]) if strand=='+' else int(p[5]); sym=p[12]
        if '_' in chrom: continue
        tss.setdefault(sym,[]).append((chrom,t))
# build per-chrom sorted TSS arrays
from collections import defaultdict
chrom_tss=defaultdict(list)
for sym,locs in tss.items():
    for c,t in locs: chrom_tss[c].append((t,sym))
for c in chrom_tss: chrom_tss[c]=sorted(chrom_tss[c])
chrom_pos={c:np.array([x[0] for x in v]) for c,v in chrom_tss.items()}
chrom_sym={c:[x[1] for x in v] for c,v in chrom_tss.items()}
print("refGene genes:",len(tss))

def genes_near(chrom,start,end,win=WIN):
    if chrom not in chrom_pos: return set()
    pos=chrom_pos[chrom]; mid=(start+end)//2
    lo=np.searchsorted(pos,mid-win); hi=np.searchsorted(pos,mid+win)
    return set(chrom_sym[chrom][lo:hi])

# ---------- TF footprint -> target genes ----------
TFS={'ESRRB_MA0141.3':'ESRRB','NANOG_NANOG':'NANOG','SOX2_MA0143.4':'SOX2','Pou5f1Sox2_MA0142.1':'OCT4-SOX2',
     'KLF4_MA0039.4':'KLF4','KLF17_MA1514.1':'KLF17','SMAD2_MA1964.1':'SMAD2/NODAL','SMAD3_MA0795.1':'SMAD3/NODAL',
     'SREBF1_MA0595.1':'SREBF1/lipid','SREBF2_MA0596.1':'SREBF2/lipid'}
tf_targets={}; tf_targets_icmgain={}
for pref,lab in TFS.items():
    fp=f'{BD}/{pref}/{pref}_overview.txt'
    if not os.path.exists(fp): continue
    d=pd.read_csv(fp,sep='\t')
    b=d[d['ICM_bound']==1]                       # bound in ICM
    tg=set()
    for _,r in b.iterrows(): tg|=genes_near(r['TFBS_chr'],r['TFBS_start'],r['TFBS_end'])
    tf_targets[lab]=tg
    # sites MORE bound in ICM than naive (footprint gained in ICM)
    bg=b[b['ICM_Naive_log2fc']>0]; tg2=set()
    for _,r in bg.iterrows(): tg2|=genes_near(r['TFBS_chr'],r['TFBS_start'],r['TFBS_end'])
    tf_targets_icmgain[lab]=tg2
    print(f"{lab:14} ICM-bound sites {len(b):6d} -> targets {len(tg):5d} | ICM-gained {len(bg):5d} -> {len(tg2):5d}")

# ---------- RNA programs: ICM vs preIm-Epi ----------
SI,C,cG,E,_,_=load_embryo_data(); ann=SI['Manual_Annotations'].astype(str).values
genes=list(C.columns)
def ln(A): A=np.asarray(A,float); lib=A.sum(1,keepdims=True); lib[lib==0]=1; return np.log1p(A/lib*1e4)
X=ln(C.values)
def program(g1,g2,topn=400):
    m1=X[ann==g1].mean(0); m2=X[ann==g2].mean(0); lfc=m1-m2
    expr=(np.maximum(m1,m2)>0.3)
    up1=np.array(genes)[expr][np.argsort(-lfc[expr])][:topn]     # up in g1
    up2=np.array(genes)[expr][np.argsort(lfc[expr])][:topn]      # up in g2
    return set(up1),set(up2)
ICM_up,EPI_up=program('ICM','preIm-Epi')
print(f"ICM-up genes {len(ICM_up)}, EPI-up genes {len(EPI_up)}")

# ---------- enrichment: TF targets vs ICM-up / EPI-up ----------
universe=set(genes)&set(tss.keys())
N=len(universe)
def enr(targets,program):
    t=targets&universe; s=program&universe; k=len(t&s)
    if len(t)==0 or len(s)==0: return 1.0,0.0,k
    p=hypergeom.sf(k-1,N,len(s),len(t))
    exp=len(t)*len(s)/N; odds=(k+.5)/(exp+.5)
    return p,odds,k
rows=[]
for lab in tf_targets:
    for prog,pset in [('ICM-up',ICM_up),('EPI-up',EPI_up)]:
        p,odds,k=enr(tf_targets[lab],pset)
        rows.append(dict(TF=lab,program=prog,p=p,odds=odds,k=k,ntarget=len(tf_targets[lab]&universe)))
R=pd.DataFrame(rows); R['neglog10p']=-np.log10(R['p']+1e-300)
print(R.pivot(index='TF',columns='program',values='odds').round(2).to_string())

# ---------- FIGURE ----------
fig,axs=plt.subplots(1,2,figsize=(15,5.8),gridspec_kw={'width_ratios':[1.15,1]})
ax=axs[0]; tfs=list(tf_targets.keys()); y=np.arange(len(tfs)); w=.38
icm_or=[R[(R.TF==t)&(R.program=='ICM-up')]['odds'].values[0] for t in tfs]
epi_or=[R[(R.TF==t)&(R.program=='EPI-up')]['odds'].values[0] for t in tfs]
icm_p =[R[(R.TF==t)&(R.program=='ICM-up')]['p'].values[0] for t in tfs]
epi_p =[R[(R.TF==t)&(R.program=='EPI-up')]['p'].values[0] for t in tfs]
order=np.argsort([np.log2(icm_or[i]/epi_or[i]) for i in range(len(tfs))])
tfs=[tfs[i] for i in order]; icm_or=[icm_or[i] for i in order]; epi_or=[epi_or[i] for i in order]
icm_p=[icm_p[i] for i in order]; epi_p=[epi_p[i] for i in order]; y=np.arange(len(tfs))
ax.barh(y+w/2,np.log2(icm_or),w,color='#c0392b',label='ICM-up program')
ax.barh(y-w/2,np.log2(epi_or),w,color='#2c3e70',label='EPI-up program (preIm epiblast)')
for i in range(len(tfs)):
    if icm_p[i]<0.05: ax.text(np.log2(icm_or[i])+.03*np.sign(np.log2(icm_or[i])+.01),i+w/2,'*',va='center',fontsize=13,weight='bold',color='#c0392b')
    if epi_p[i]<0.05: ax.text(np.log2(epi_or[i])+.03*np.sign(np.log2(epi_or[i])+.01),i-w/2,'*',va='center',fontsize=13,weight='bold',color='#2c3e70')
ax.set_yticks(y); ax.set_yticklabels(tfs,fontsize=10); ax.axvline(0,color='k',lw=.6)
ax.set_xlabel('footprint-target enrichment  log2(odds ratio)',fontsize=10.5); ax.legend(fontsize=9,loc='lower right')
ax.set_title('A   Which TF footprint-targets drive each side of the ICM→epiblast transition\n(* hypergeometric p<0.05)',fontsize=11,loc='left',weight='bold')
# B: diverging net ICM- vs EPI-driver preference = log2(odds_ICM / odds_EPI)
ax=axs[1]
net=[(t,np.log2(R[(R.TF==t)&(R.program=='ICM-up')]['odds'].values[0]/R[(R.TF==t)&(R.program=='EPI-up')]['odds'].values[0])) for t in tf_targets]
net=sorted(net,key=lambda x:x[1]); nl=[x[0] for x in net]; nv=[x[1] for x in net]
colb=['#c0392b' if v>0 else '#2c3e70' for v in nv]
ax.barh(range(len(nl)),nv,color=colb,edgecolor='k',lw=.4)
ax.axvline(0,color='k',lw=.8)
ax.set_yticks(range(len(nl))); ax.set_yticklabels(nl,fontsize=10)
ax.set_xlabel('net preference  log2(ICM-up odds / EPI-up odds)',fontsize=10)
ax.text(.98,.06,'← EPI-driver',transform=ax.transAxes,ha='right',fontsize=9,color='#2c3e70',weight='bold')
ax.text(.98,.12,'ICM-driver →',transform=ax.transAxes,ha='right',fontsize=9,color='#c0392b',weight='bold')
ax.set_title('B   Net readout: OCT4-SOX2 / NANOG / ESRRB are ICM-drivers;\nNODAL/SMAD & lipid (SREBF) are epiblast-drivers',fontsize=11,loc='left',weight='bold')
fig.suptitle('Footprint-anchored regulatory readout: TOBIAS TF targets vs the ICM↔epiblast RNA program (GSE101571 ATAC, 3PN caveat)',fontsize=11.5,y=1.02)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_GroundState_FootprintNetwork.png',dpi=175,bbox_inches='tight')
print("SAVED Fig_GroundState_FootprintNetwork.png")
R.to_csv('/tmp/footprint_enrichment.csv',index=False)
json.dump({lab:{'targets':len(tf_targets[lab]),'icm_up_p':float(R[(R.TF==lab)&(R.program=='ICM-up')]['p'].values[0]),
                'epi_up_p':float(R[(R.TF==lab)&(R.program=='EPI-up')]['p'].values[0])} for lab in tf_targets},
          open('/tmp/footprint_network.json','w'),indent=2)
print("DONE footprint network")
