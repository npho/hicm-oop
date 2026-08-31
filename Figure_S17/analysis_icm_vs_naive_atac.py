# ============================================================================
# Figure S17  --  ICM vs naive ATAC comparison
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   ICM versus naive-hESC (epiblast proxy) accessibility.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   ICM vs naive ATAC comparison
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S17. ICM-versus-naïve-hESC (epiblast proxy) accessibility. The
#   confident, promoter-anchored normal-ploidy ICM peaks are shared with the
#   naïve/epiblast-like state, consistent with a shared (not ICM-specific)
#   circuit. No human embryo epiblast ATAC exists; naïve hESC is used as the
#   epiblast proxy.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_icm_vs_naive_atac.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
(b) ICM-vs-naive-hESC accessibility proxy: are normal-ploidy (2PN) ICM-accessible regions
SHARED with the naive/epiblast-like state or ICM-specific? Naive hESC = in-vitro epiblast proxy
(no human embryo epiblast ATAC exists). Classify each 2PN ICM peak by naive-hESC signal vs a
random background; report shared vs ICM-specific fraction and where ground-state loci fall.
"""
import os,sys,gzip,json,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyBigWig
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE='/Users/adelosan/cESFW_analysis'; os.chdir(HERE); rng=np.random.RandomState(0)
PLOTS=f'{HERE}/plots/revision_results'; STD=['chr%d'%i for i in range(1,23)]+['chrX']
bwN=pyBigWig.open('atac_rt_pipeline/tobias_out/corr_Naive/Naive_uncorrected.bw')
bwI=pyBigWig.open('atac_rt_pipeline/tobias_out/corr_ICM/ICM_uncorrected.bw')

peaks=[]
with gzip.open('data/gao2018_atac/icm_2pn_peaks.bed.gz','rt') as f:
    for line in f:
        p=line.split('\t')
        if p[0] in STD: mid=(int(p[1])+int(p[2]))//2; peaks.append((p[0],mid-250,mid+250))
P=pd.DataFrame(peaks,columns=['chr','start','end'])
def shuffle(df):
    out=[]
    for c,s,e in df[['chr','start','end']].values:
        if c in bwN.chroms(): st=rng.randint(1,bwN.chroms()[c]-600); out.append((c,st,st+(e-s)))
    return pd.DataFrame(out,columns=['chr','start','end'])
ctrl=shuffle(P)
def sig(bw,df):
    v=np.zeros(len(df))
    for i,(c,s,e) in enumerate(df[['chr','start','end']].values):
        if c in bw.chroms() and 0<=s and e<=bw.chroms()[c]:
            m=bw.stats(c,int(s),int(e),type='mean')[0]
            if m is not None: v[i]=m
    return v
nP=sig(bwN,P); nC=sig(bwN,ctrl)
thr=np.percentile(nC,95)                 # naive-accessible if above 95th pct of random naive signal
shared=nP>thr
print(f"2PN ICM peaks {len(P)} | naive-accessible (shared, epiblast-like) {shared.mean()*100:.0f}% | ICM-specific {(~shared).mean()*100:.0f}% (naive 95pct thr={thr:.2f})")

# assign shared vs ICM-specific peaks to genes; check ground-state/conserved loci
def load_tss(fp):
    t=defaultdict(list)
    for l in gzip.open(fp,'rt'):
        p=l.rstrip().split('\t'); c=p[2]
        if '_' in c: continue
        t[p[12]].append((c,int(p[4]) if p[3]=='+' else int(p[5])))
    return t
ht=load_tss('data/GSE101571/refGene_hg19.txt.gz')
gene_class={}
for i in range(len(P)):
    c,s,e=P.iloc[i][['chr','start','end']]; mid=(s+e)//2
GS=['ESRRB','TFCP2L1','KLF2','KLF4','KLF5','TBX3','PRDM14','ZFP42','NR0B1','DPPA3','KLF17','SALL4','NANOG','POU5F1','SOX2','GATA6','GDF3','TDGF1','FGF4']
# for each GS gene, is its TSS near a shared or ICM-specific ICM peak?
tss_g={g:ht[g][0] for g in GS if g in ht}
res=[]
for g,(c,t) in tss_g.items():
    near=[(shared[i]) for i in range(len(P)) if P.iloc[i]['chr']==c and abs((P.iloc[i]['start']+P.iloc[i]['end'])//2 - t)<20000]
    if near: res.append((g,'shared(naive+)' if any(near) else 'ICM-specific','; peaks='+str(len(near))))
    else: res.append((g,'no 2PN ICM peak',''))
for g,cl,x in res: print(f"  {g:8} {cl}")

# ---- figure ----
fig,axs=plt.subplots(1,2,figsize=(12.5,5),gridspec_kw={'width_ratios':[1,1.1]})
ax=axs[0]
ax.hist(np.log1p(nP),bins=60,alpha=.7,color='#c0392b',label='naive signal @ ICM peaks',density=True)
ax.hist(np.log1p(nC),bins=60,alpha=.5,color='#999',label='naive signal @ random',density=True)
ax.axvline(np.log1p(thr),color='k',ls='--',lw=1); ax.set_xlabel('log(1+ naive-hESC accessibility)',fontsize=10)
ax.set_ylabel('density',fontsize=10); ax.legend(fontsize=8.5)
ax.set_title(f'A   {shared.mean()*100:.0f}% of 2PN ICM peaks are naive-accessible (shared)',fontsize=11,loc='left',weight='bold')
ax=axs[1]; lab=['shared\n(naive+ / epiblast-like)','ICM-specific\n(naive-low)']
ax.bar([0,1],[shared.mean()*100,(~shared).mean()*100],color=['#2c3e70','#c0392b'])
for i,v in enumerate([shared.mean()*100,(~shared).mean()*100]): ax.text(i,v+1,f'{v:.0f}%',ha='center',fontsize=11,weight='bold')
ax.set_xticks([0,1]); ax.set_xticklabels(lab,fontsize=9.5); ax.set_ylabel('% of 2PN ICM peaks',fontsize=10); ax.set_ylim(0,100)
ax.set_title('B   ICM chromatin is largely shared with the naive/epiblast-like\nstate — consistent with the scRNA (circuit is shared, not ICM-specific)',fontsize=10,loc='left',weight='bold')
fig.suptitle('ICM vs naive-hESC (epiblast proxy) accessibility — no human embryo epiblast ATAC exists',fontsize=11.5,y=1.02)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_ICM_vs_Naive_ATAC.png',dpi=170,bbox_inches='tight')
print("SAVED Fig_ICM_vs_Naive_ATAC.png")
json.dump({'shared_pct':float(shared.mean()*100),'icm_specific_pct':float((~shared).mean()*100),'n_peaks':len(P)},
          open('/tmp/icm_vs_naive.json','w'),indent=2)
print("DONE")
