# ============================================================================
# Figure S19  --  per-study stage-sequence replication
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Inter-study robustness of every key ICM analysis, ordered by main
#   figure.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   per-study stage-sequence replication
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S19. Inter-study robustness of the key ICM analyses, ordered by
#   main figure (Petropoulos, Yanagida, Meistermann; Yanagida-inclusive
#   atlas; epiblast = preIm-Epi and Epi/Hypo-Branch cells; per-study n shown
#   where relevant). (A) Overview: ICM→epiblast expression change per study;
#   identical sign across columns indicates a pooling-robust trend. (B, C)
#   Figure 2: (B) ICM methylation machinery, UHRF1 near-absent in all
#   studies (demethylation-permissive), TET2 present; (C) c-MYC low across
#   preimplantation (ICM ≈ epiblast) and rising only post-implantation.
#   (D–G) Figure 3: (D) ESRRB declines to the epiblast in every study while
#   TFCP2L1 is switched on at the ICM and maintained; (E) the
#   OCT4/SOX2/SALL4 core is maintained (Δ≈0) while NANOG rises; (F) the
#   circuit is pluripotency-directed (pluripotency − TE > 0;
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/stage_anatomy/epi_replicate.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
import decoupler as dc

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
Counts=Counts.fillna(0); allg=list(Counts.columns)
STAGES=["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi"]; LAB=["ICM/TE","ICM","Epi/Hypo Br","EPI"]
net=pd.read_csv("/Users/adelosan/cesfw_runs/collectri_net.csv")
gmt={}
for line in open("/Users/adelosan/.cache/gseapy/Enrichr.MSigDB_Hallmark_2020.gmt"):
    pp=line.strip().split("\t"); gmt[pp[0]]=[g for g in pp[2:] if g]

def analyze(ds):
    d=Sample_Info[(Sample_Info.Dataset==ds)&(Sample_Info.Manual_Annotations.isin(STAGES))]
    ncells={s:(d.Manual_Annotations==s).sum() for s in STAGES}
    usable=[s for s in STAGES if ncells[s]>=3]
    if len(usable)<3: return None,ncells
    cs=Counts.loc[d.index]; lib=cs.sum(axis=1); lib[lib==0]=1; lc=np.log2(cs.div(lib,axis=0)*1e6+1)
    stage=d["Manual_Annotations"].values
    def delta(genes):
        p=[g for g in genes if g in lc.columns]; bg=[g for g in allg if g not in p]
        return lc[p].mean(axis=1)-lc[bg].mean(axis=1)
    chol=[delta(gmt["Cholesterol Homeostasis"])[stage==s].mean() if ncells[s]>=3 else np.nan for s in STAGES]
    acts=dc.mt.ulm(data=lc,net=net,tmin=5); tf=acts[0] if isinstance(acts,tuple) else acts
    def tfstage(name):
        return [tf.loc[d.index[stage==s],name].mean() if (ncells[s]>=3 and name in tf.columns) else np.nan for s in STAGES]
    return {"chol":chol,"SREBF1":tfstage("SREBF1"),"SREBF2":tfstage("SREBF2"),
            "FOXO3":tfstage("FOXO3"),"POU5F1":tfstage("POU5F1"),"MYC":tfstage("MYC")}, ncells

print("Replication of the EPI-trajectory sequence: values per stage [ICM/TE, ICM, Epi/Hypo Br, EPI]")
print("Key claims: (1) cholesterol/SREBF RISE across the trajectory; (2) SREBF1/2 jump is largest")
print("            entering Epi/Hypo Branch; (3) POU5F1 highest at EPI; (4) MYC high early.\n")
for ds in ["Petropoulos 2016","Yanagida 2020","Meistermann 2021"]:
    r,nc=analyze(ds)
    print(f"### {ds}  (n/stage: {nc}) ###")
    if r is None: print("  too few cells in >1 stage\n"); continue
    for k in ["chol","SREBF1","SREBF2","FOXO3","POU5F1","MYC"]:
        v=r[k]; vs=" ".join(f"{x:+.2f}" if not np.isnan(x) else "  NA " for x in v)
        # largest step-up location
        steps=[v[i+1]-v[i] if not(np.isnan(v[i])or np.isnan(v[i+1])) else -np.inf for i in range(3)]
        loc=["->ICM","->EpiHypoBr","->EPI"][int(np.argmax(steps))] if max(steps)>-np.inf else "NA"
        print(f"  {k:8s}: {vs}   | biggest rise at {loc}")
    print()
print("DONE")
