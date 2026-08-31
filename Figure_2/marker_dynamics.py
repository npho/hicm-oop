# ============================================================================
# Figure 2  --  naive-marker dynamics across stages
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The canonical naive network is assembled only AFTER the ICM: ICM is
#   demethylation-permissive (high TET2, near-absent UHRF1) and low MYC.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   naive-marker dynamics across stages
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/stage_anatomy/marker_dynamics.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
import warnings; warnings.filterwarnings("ignore")
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data(); Counts=Counts.fillna(0); allg=list(Counts.columns)
STAGES=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi"]; LAB=["ICM/TE","ICM","Epi/Hypo\nBr","EPI"]
gmt={}
for line in open("/Users/adelosan/.cache/gseapy/Enrichr.MSigDB_Hallmark_2020.gmt"):
    pp=line.strip().split("\t"); gmt[pp[0]]=[g for g in pp[2:] if g]
DS=[("Petropoulos 2016","#16a34a"),("Yanagida 2020","#3b7dd8"),("Meistermann 2021","#f59e0b")]
def stage_vals(ds,gene_or_path,is_path=False):
    d=Sample_Info[(Sample_Info.Dataset==ds)&(Sample_Info.Manual_Annotations.isin(STAGES))]
    cs=Counts.loc[d.index]; lib=cs.sum(axis=1); lib[lib==0]=1; lc=np.log2(cs.div(lib,axis=0)*1e6+1)
    stage=d["Manual_Annotations"].values; nc={s:(stage==s).sum() for s in STAGES}
    if is_path:
        p=[x for x in gmt[gene_or_path] if x in lc.columns]; bg=[x for x in allg if x not in p]
        sc=lc[p].mean(axis=1)-lc[bg].mean(axis=1)
        return [sc[stage==s].mean() if nc[s]>=3 else np.nan for s in STAGES]
    if gene_or_path not in lc.columns: return [np.nan]*4
    return [lc.loc[d.index[stage==s],gene_or_path].mean() if nc[s]>=3 else np.nan for s in STAGES]

panels=[("ESRRB","ESRRB (pre-naive ICM marker)",False),("NANOG","NANOG (naive epiblast marker)",False),
        ("KLF17","KLF17 (biphasic)",False),("SOX2","SOX2",False),
        ("NODAL","NODAL (TGF-β ligand)",False),("TGF-beta Signaling","TGF-β signalling (Hallmark)",True)]
fig,axes=plt.subplots(2,3,figsize=(14,8))
for ax,(g,title,isp) in zip(axes.flat,panels):
    for ds,c in DS:
        v=stage_vals(ds,g,isp)
        ax.plot(range(4),v,marker="o",color=c,lw=2,label=ds.split()[0])
    ax.set_xticks(range(4)); ax.set_xticklabels(LAB,fontsize=8)
    ax.set_title(title,fontsize=10,fontweight="bold")
    ax.set_ylabel("delta-score" if isp else "mean log2CPM",fontsize=8)
axes.flat[0].legend(fontsize=7)
plt.suptitle("Marker and signalling dynamics along the epiblast trajectory (three datasets)\nESRRB marks the PRE-NAIVE ICM; the naive network (NANOG/SOX2/NODAL/TGF-β) switches on at the naive epiblast",
             fontweight="bold",fontsize=12)
plt.tight_layout(); plt.savefig("/Users/adelosan/cesfw_runs/manuscript_build/figures/fig15_markers.png",dpi=150,bbox_inches="tight")
print("saved fig15_markers.png")
