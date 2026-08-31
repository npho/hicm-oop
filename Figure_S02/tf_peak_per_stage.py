# ============================================================================
# Figure S2  --  TF activity peak per stage
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Candidate stage-specific regulators along the epiblast trajectory.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   TF activity peak per stage
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S2. Candidate stage-specific regulators along the epiblast
#   trajectory. Stage-specific transcription-factor activity (row z-scored;
#   top replicated candidates). The ICM/TE Branch is enriched for the naïve
#   factor KLF17, the ICM for the metabolic regulator PPARD, and the
#   preimplantation epiblast for NOTCH1, HEY2 and the Wnt effector TCF4.
#   These are candidates for future functional testing, not validated
#   drivers.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/stage_anatomy/tf_peak_per_stage.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
import decoupler as dc
import warnings; warnings.filterwarnings("ignore")

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data(); Counts=Counts.fillna(0)
STAGES=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi"]; LAB=["ICM/TE","ICM","Epi/Hypo Br","EPI"]
net=pd.read_csv("/Users/adelosan/cesfw_runs/collectri_net.csv")

def tf_by_stage(ds):
    d=Sample_Info[(Sample_Info.Dataset==ds)&(Sample_Info.Manual_Annotations.isin(STAGES))]
    cs=Counts.loc[d.index]; lib=cs.sum(axis=1); lib[lib==0]=1; lc=np.log2(cs.div(lib,axis=0)*1e6+1)
    stage=d["Manual_Annotations"].values
    acts=dc.mt.ulm(data=lc,net=net,tmin=5); tf=acts[0] if isinstance(acts,tuple) else acts
    nc={s:(stage==s).sum() for s in STAGES}
    m=pd.DataFrame({LAB[i]:(tf[stage==STAGES[i]].mean(axis=0) if nc[STAGES[i]]>=3 else np.nan) for i in range(4)})
    return m,nc

petro,nc=tf_by_stage("Petropoulos 2016")
print("Petropoulos cells/stage:",nc)
# z-score across stages so 'peak' means stage-specific, not just globally high
z=petro.sub(petro.mean(axis=1),axis=0).div(petro.std(axis=1)+1e-9,axis=0)
petro["peak"]=petro.idxmax(axis=1); petro["peak_z"]=z.max(axis=1); petro["spread"]=petro[LAB].max(axis=1)-petro[LAB].min(axis=1)

# candidate regulators = TFs that peak at each stage, ranked by stage-specificity (peak z) among reasonably-active TFs
print("\n"+"="*70); print("CANDIDATE STAGE-SPECIFIC REGULATORS (TFs peaking at each stage)"); print("="*70)
per_stage={}
for l in LAB:
    sub=petro[(petro["peak"]==l)&(petro["spread"]>0.3)].sort_values("peak_z",ascending=False)
    per_stage[l]=sub.head(12).index.tolist()
    print(f"\n### {l} ###")
    print("  "+", ".join(f"{g}" for g in sub.head(12).index))

# replication: for each candidate, does it also peak at same stage in Yanagida?
yana,ncy=tf_by_stage("Yanagida 2020")
yana_peak=yana.idxmax(axis=1)
print("\n"+"="*70); print("REPLICATION in Yanagida (does candidate peak at same stage?)"); print("="*70)
for l in LAB:
    cands=per_stage[l]
    rep=[g for g in cands if g in yana_peak.index and yana_peak[g]==l]
    print(f"  {l:12s}: {len(rep)}/{len(cands)} replicate -> {', '.join(rep) if rep else '(none)'}")

petro[LAB+["peak","peak_z","spread"]].to_csv("/Users/adelosan/cesfw_runs/stage_anatomy/tf_peak_per_stage.csv")

# figure: top candidates per stage
fig,ax=plt.subplots(figsize=(7,9))
order=[]; labels=[]
for l in LAB:
    for g in per_stage[l][:8]: order.append(g); labels.append(l)
zz=z.loc[order]
im=ax.imshow(zz.values,aspect="auto",cmap="RdBu_r",vmin=-1.4,vmax=1.4)
ax.set_xticks(range(4)); ax.set_xticklabels(LAB,rotation=30,ha="right")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order,fontsize=7)
# stage separators
cum=0
for l in LAB:
    n=len(per_stage[l][:8]); ax.axhline(cum-0.5,color="black",lw=1); cum+=n
plt.colorbar(im,ax=ax,shrink=0.3,label="row z")
ax.set_title("Candidate stage-specific regulators (TFs peaking at each stage)")
plt.tight_layout(); plt.savefig("/Users/adelosan/cesfw_runs/stage_anatomy/tf_peak_per_stage.png",dpi=150,bbox_inches="tight")
print("\nSaved tf_peak_per_stage.png"); print("DONE")
