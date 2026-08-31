# ============================================================================
# Figure 3  --  G,H,I: normal-ploidy ICM chromatin accessibility
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM ignites a conserved ground-state circuit (ESRRB/TFCP2L1) on a
#   maintained OCT4/SOX2/NANOG core, with mouse ChIP conservation and ICM
#   footprints.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   G,H,I: normal-ploidy ICM chromatin accessibility
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_icm_chromatin_priming.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
FIX the 3PN caveat + build the naive-vs-primed chromatin test (proposal Fig 4A), using
NORMAL-PLOIDY (2PN) human ICM peaks (Gao 2018, GSE101571) to define regions and hESC
accessibility (naive/primed, ploidy-irrelevant) quantified there.
 A  ploidy concordance: 2PN vs 3PN ICM peaks + 3PN signal elevated at 2PN peaks (not a ploidy artifact)
 B  naive vs primed accessibility enrichment at 2PN ICM peaks (each vs its own random background)
 C  ground-state loci: naive vs primed signal at their 2PN-ICM peaks
"""
import os,sys,gzip,json,warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyBigWig
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
HERE='/Users/adelosan/cESFW_analysis'; os.chdir(HERE); rng=np.random.RandomState(0)
PLOTS=f'{HERE}/plots/revision_results'
GAO='data/gao2018_atac'; TB='atac_rt_pipeline/tobias_out'
STDCHR=['chr%d'%i for i in range(1,23)]+['chrX']

def load_bed(fp):
    rows=[]
    with gzip.open(fp,'rt') as f:
        for line in f:
            p=line.split('\t')
            if p[0] in STDCHR: rows.append((p[0],int(p[1]),int(p[2])))
    return pd.DataFrame(rows,columns=['chr','start','end'])
icm2=load_bed(f'{GAO}/icm_2pn_peaks.bed.gz'); icm3=load_bed(f'{GAO}/icm_3pn_peaks.bed.gz')
print(f"2PN ICM peaks {len(icm2)} | 3PN ICM peaks {len(icm3)}")

# ---- ploidy concordance: fraction of 3PN peaks overlapping a 2PN peak ----
def overlap_frac(A,B):
    hit=0
    for c in STDCHR:
        a=A[A['chr']==c][['start','end']].values; b=B[B['chr']==c][['start','end']].values
        if len(a)==0 or len(b)==0: continue
        bs=b[:,0]; be=b[:,1]; order=np.argsort(bs); bs=bs[order]; be=be[order]
        for s,e in a:
            i=np.searchsorted(bs,e)
            # any b with bs<e and be>s among [max(0,i-50):i]
            lo=max(0,i-200)
            if np.any((bs[lo:i]<e)&(be[lo:i]>s)): hit+=1
    return hit/len(A)
f_3in2=overlap_frac(icm3,icm2); f_2in3=overlap_frac(icm2,icm3)
print(f"3PN peaks recovered in 2PN: {100*f_3in2:.1f}% | 2PN peaks in 3PN: {100*f_2in3:.1f}%")

# ---- control regions: shuffle 2PN peaks within chrom ----
bwN=pyBigWig.open(f'{TB}/corr_Naive/Naive_uncorrected.bw')
bwI=pyBigWig.open(f'{TB}/corr_ICM/ICM_uncorrected.bw')
chroms=bwN.chroms()
def shuffle(df):
    out=[]
    for _,r in df.iterrows():
        L=r['end']-r['start']; c=r['chr']
        if c not in chroms: continue
        s=rng.randint(1,max(2,chroms[c]-L-1)); out.append((c,s,s+L))
    return pd.DataFrame(out,columns=['chr','start','end'])
def summit_win(df,w=250):
    mid=((df['start']+df['end'])//2).values
    return pd.DataFrame({'chr':df['chr'].values,'start':np.maximum(1,mid-w),'end':mid+w})
icm2=summit_win(icm2)          # summit-centered 500bp windows (broad-peak dilution fix)
ctrl=shuffle(icm2)

def bw_mean(bw,df):             # None/no-coverage -> 0.0 (sparse bigwig; NaN-drop was the bug)
    v=np.zeros(len(df))
    for i,(c,s,e) in enumerate(df[['chr','start','end']].values):
        if c in bw.chroms() and e<=bw.chroms()[c]:
            m=bw.stats(c,int(s),int(e),type='mean')[0]
            if m is not None: v[i]=m
    return v
naive_icm=bw_mean(bwN,icm2); naive_ctrl=bw_mean(bwN,ctrl)
icm_icm  =bw_mean(bwI,icm2); icm_ctrl =bw_mean(bwI,ctrl)

# ---- primed signal from wig (parse once, assign to icm2 + ctrl) ----
def wig_means(regions):
    # per-chrom sorted region index
    idx={};
    for c in STDCHR:
        r=regions[regions['chr']==c][['start','end']].values
        order=np.argsort(r[:,0]) if len(r) else np.array([],int)
        idx[c]=(r[order] if len(r) else r, np.where(regions['chr'].values==c)[0][order] if len(r) else np.array([],int))
    acc=np.zeros(len(regions)); cnt=np.zeros(len(regions))
    with gzip.open('data/GSE101571/primed_hESC_rpkm.wig.gz','rt') as f:
        for line in f:
            p=line.split('\t')
            c=p[0]
            if c not in idx or len(idx[c][0])==0: continue
            s=int(p[1]); val=float(p[3]); starts=idx[c][0][:,0]; ends=idx[c][0][:,1]; gi=idx[c][1]
            j=np.searchsorted(starts,s+100)
            lo=max(0,j-50)
            for k in range(lo,min(j+1,len(starts))):
                if starts[k]<s+100 and ends[k]>s:
                    acc[gi[k]]+=val; cnt[gi[k]]+=1
    return np.where(cnt>0,acc/np.maximum(cnt,1),0.0)   # no-coverage -> 0 (comparable to bigwig zero-fill)
print("parsing primed wig (both region sets)...")
allreg=pd.concat([icm2.assign(g='icm'),ctrl.assign(g='ctrl')],ignore_index=True)
pv=wig_means(allreg)
primed_icm=pv[:len(icm2)]; primed_ctrl=pv[len(icm2):]

def enr(peak,ctrl):
    p=np.nanmean(peak); c=np.nanmean(ctrl); return p,c,p/c if c>0 else np.nan
en_naive=enr(naive_icm,naive_ctrl); en_primed=enr(primed_icm,primed_ctrl); en_icm=enr(icm_icm,icm_ctrl)
print(f"enrichment at 2PN ICM peaks vs random: naive {en_naive[2]:.2f}x | primed {en_primed[2]:.2f}x | ICM(3PN) {en_icm[2]:.2f}x")

# ---- ground-state loci: TSS from refGene, nearest 2PN peak + naive/primed signal ----
tss={}
with gzip.open('data/GSE101571/refGene_hg19.txt.gz','rt') as f:
    for line in f:
        p=line.rstrip().split('\t'); c=p[2]; st=p[3]; t=int(p[4]) if st=='+' else int(p[5]); tss.setdefault(p[12],[]).append((c,t))
GS=['ESRRB','TFCP2L1','NANOG','KLF17','KLF4','NODAL','TDGF1','LEFTY1','GATA6','POU5F1','SOX2']
loci=[]
for g in GS:
    if g not in tss: continue
    c,t=tss[g][0]
    pk=icm2[(icm2['chr']==c)&(icm2['start']<t+20000)&(icm2['end']>t-20000)]
    has=len(pk)>0
    reg=pd.DataFrame([[c,max(1,t-1000),t+1000]],columns=['chr','start','end'])
    nv=bw_mean(bwN,reg)[0]; pr=wig_means(reg)[0]
    loci.append(dict(gene=g,has_2PN_peak=has,naive=nv,primed=pr))
L=pd.DataFrame(loci); print(L.to_string())

# ---- FIGURE ----
fig,axs=plt.subplots(1,3,figsize=(16.5,5.2),gridspec_kw={'width_ratios':[0.8,1,1.25]})
# A ploidy concordance
ax=axs[0]; ax.bar([0,1],[100*f_3in2,100*f_2in3],color=['#c0392b','#e67e22'])
ax.set_xticks([0,1]); ax.set_xticklabels(['3PN peaks\nrecovered in 2PN','2PN peaks\nin 3PN'],fontsize=9)
ax.set_ylabel('% overlap',fontsize=11); ax.set_ylim(0,100)
for i,v in enumerate([f_3in2,f_2in3]): ax.text(i,100*v+2,f'{100*v:.0f}%',ha='center',fontsize=10,weight='bold')
ax.set_title('G   Normal-ploidy (2PN) and 3PN ICM\npeaks concur → not a ploidy artifact',fontsize=10.5,loc='left',weight='bold')
# B enrichment
ax=axs[1]; vals=[en_naive[2],en_primed[2],en_icm[2]]; cols=['#2c3e70','#7f8c8d','#c0392b']
ax.bar(range(3),vals,color=cols); ax.axhline(1,color='k',ls='--',lw=1)
ax.set_xticks(range(3)); ax.set_xticklabels(['naïve\nhESC','primed\nhESC','ICM (3PN)\nself-check'],fontsize=9)
ax.set_ylabel('accessibility enrichment\nat 2PN ICM peaks (vs random)',fontsize=10.5)
for i,v in enumerate(vals): ax.text(i,v+.05*max(vals),f'{v:.1f}x',ha='center',fontsize=10,weight='bold')
ax.set_title('H   Normal-ploidy ICM-accessible regions are\nnaïve-accessible >> primed-accessible',fontsize=10.5,loc='left',weight='bold')
# C loci
ax=axs[2]; gg=L['gene'].values; y=np.arange(len(gg)); w=.38
nn=L['naive'].values/np.nanmax(L['naive'].values); pp=L['primed'].values/np.nanmax(L['primed'].values)
ax.barh(y+w/2,nn,w,color='#2c3e70',label='naïve hESC'); ax.barh(y-w/2,pp,w,color='#7f8c8d',label='primed hESC')
for i,g in enumerate(gg):
    if L['has_2PN_peak'].values[i]: ax.text(-0.02,i,'●',ha='right',va='center',color='#c0392b',fontsize=9)
ax.set_yticks(y); ax.set_yticklabels(gg,fontsize=9); ax.set_xlabel('normalized accessibility at TSS±1kb',fontsize=10)
ax.legend(fontsize=8,loc='lower right'); ax.set_title('I   Ground-state loci (● = has normal-ploidy 2PN ICM peak):\nnaïve-accessible in hESC',fontsize=10,loc='left',weight='bold')
fig.suptitle('ICM chromatin is primed for naïve pluripotency — on NORMAL-PLOIDY (2PN) ICM peaks (Gao 2018, GSE101571)',fontsize=12.5,y=1.02)
plt.tight_layout(); fig.savefig(f'{PLOTS}/Fig_ICM_ChromatinPriming.png',dpi=170,bbox_inches='tight')
print("SAVED Fig_ICM_ChromatinPriming.png")
json.dump({'peaks_2PN':len(icm2),'peaks_3PN':len(icm3),'concord_3in2':f_3in2,'concord_2in3':f_2in3,
           'enrich_naive':en_naive[2],'enrich_primed':en_primed[2],'enrich_ICM3PN':en_icm[2],
           'loci':L.to_dict('records')},open('/tmp/icm_chromatin_priming.json','w'),default=float,indent=2)
print("DONE")
