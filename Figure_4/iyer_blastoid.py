# ============================================================================
# Figure 4  --  Iyer blastoid Ctrl vs RapaLink (mTOR inhibition does not induce ESRRB)
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   No in-vitro system stabilises the ICM programme, but rare genuine ESRRB+
#   cells arise. Two independent measurements (reference mapping; direct
#   signature scoring).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   Iyer blastoid Ctrl vs RapaLink (mTOR inhibition does not induce ESRRB)
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure 4. No in-vitro system builds the ICM-like state, but rare genuine
#   ESRRB+ cells arise. (A, B) cESFW reference-mapping of naive/reset PSCs
#   and blastoids onto the human embryo atlas (cosine k-NN on cESFW genes;
#   leave-one-out ICM recall 55%); Naive PSCs map predominantly to the
#   preimplantation epiblast (or, for e4CL/8CLC, the 8-cell/morula stage),
#   whereas blastoids map to epiblast and trophectoderm. The embryo profiled
#   alongside serves as a positive control, with 22.8% ICM. (C) Reference-
#   mapped (Measurement 1): all naive-PSC and blastoid systems contain 2.5%
#   or less ICM versus 22.8% for the embryo. (D) Direct scoring of the dual-
#   competent ICM signature (Measurement 2): fraction of cells co-expressing
#   ESRRB+GATA6+ (solid) and ESRRB+TFCP2L1+GATA6+ (hatched). 93% of cESFW
#   embryo ICM-family cells versus under 5% in every steady-state in-vitro
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/iyer_blastoid.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os, gzip
from scipy.io import mmread
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import mannwhitneyu

OUT="/Users/adelosan/cesfw_runs/iyer_blastoid/"; os.makedirs(OUT,exist_ok=True)
IY="/Users/adelosan/cESFW_analysis/data/Iyer2024_Blastoid"
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data(); Counts=Counts.fillna(0); allg=list(Counts.columns)
feat=[g for g in cESFW_Genes if g in Counts.columns]
STAGES=["Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
def l2cpm(cs):
    lib=cs.sum(axis=1); lib[lib==0]=1; return np.log2(cs.div(lib,axis=0)*1e6+1)
ref=Sample_Info.index[Sample_Info["Manual_Annotations"].isin(STAGES)]
ref_lc=l2cpm(Counts.loc[ref])[feat].fillna(0)
ref_lab=Sample_Info.loc[ref,"Manual_Annotations"].values
KNN=KNeighborsClassifier(n_neighbors=20,metric="correlation",weights="distance").fit(ref_lc.values,ref_lab)

def load_10x(prefix):
    m=mmread(f"{IY}/{prefix}_matrix.mtx.gz").tocsc()  # genes x cells
    feats=pd.read_csv(f"{IY}/{prefix}_features.tsv.gz",sep="\t",header=None)
    genes=feats[1].values if feats.shape[1]>1 else feats[0].values
    bc=pd.read_csv(f"{IY}/{prefix}_barcodes.tsv.gz",header=None)[0].values
    df=pd.DataFrame.sparse.from_spmatrix(m.T, index=bc, columns=genes)  # cells x genes
    return df
meta=pd.read_csv(f"{IY}/preimpl_meta.tsv.gz",sep="\t")
print("blastoid celltypes:",meta.celltype.value_counts().to_dict())

gmt={}
for line in open("/Users/adelosan/.cache/gseapy/Enrichr.MSigDB_Hallmark_2020.gmt"):
    pp=line.strip().split("\t"); gmt[pp[0]]=[g for g in pp[2:] if g]
srebf=["HMGCR","HMGCS1","SQLE","LDLR","FDFT1","DHCR7","MVD","MVK","FASN","ACACA","SCD"]
def delta(lc,genes):
    p=[g for g in genes if g in lc.columns]; bg=[g for g in allg if g not in p and g in lc.columns]
    return lc[p].mean(axis=1)-lc[bg].mean(axis=1)

results={}
for prefix,cond in [("ctrl","Ctrl"),("rapa","RapaLink")]:
    df=load_10x(prefix)
    df=df.loc[:, ~df.columns.duplicated()]
    # match to meta (EPI cells of this condition)
    mm=meta[meta.CellID.str.contains("267607" if cond=="Ctrl" else "267608")]
    # map barcodes: meta CellID like 'Ctrl_96h-267607_AATC...-1'; matrix barcodes like 'AATC...-1'
    mm=mm.copy(); mm["bc"]=mm.CellID.str.split("_").str[-1]
    common=[b for b in mm.bc if b in df.index]
    sub=df.loc[common]; ct=mm.set_index("bc").loc[common,"celltype"].values
    dcpm=df.astype(float)
    lib=sub.astype(float).sum(axis=1); lib[lib==0]=1
    lc=np.log2(sub.astype(float).div(lib,axis=0)*1e6+1)
    # project
    X=pd.DataFrame(0.0,index=lc.index,columns=feat); pres=[g for g in feat if g in lc.columns]
    X[pres]=lc[pres].fillna(0).values
    pred=KNN.predict(X.values)
    results[cond]={"pred":pred,"ct":ct,"lc":lc,
                   "ESRRB":lc["ESRRB"].values if "ESRRB" in lc.columns else np.full(len(lc),np.nan),
                   "mTOR":delta(lc,gmt["mTORC1 Signaling"]).values,
                   "sterol":delta(lc,srebf).values}
    epi=ct=="EPI"
    vc=pd.Series(pred[epi]).value_counts(normalize=True)
    print(f"\n### {cond} blastoid EPI cells (n={epi.sum()}) — embryo-stage mapping ###")
    print("  "+", ".join(f"{s}:{100*vc.get(s,0):.0f}%" for s in STAGES if vc.get(s,0)>=0.03))
    print(f"  % mapping to pre-naive ICM (ICM+ICM/TE): {100*np.isin(pred[epi],['ICM','ICM/TE Branch']).mean():.0f}%")

# ctrl vs rapa: ESRRB, mTOR, sterol in EPI cells
print("\n"+"="*66); print("mTOR PERTURBATION TEST (RapaLink vs Ctrl, blastoid EPI cells)"); print("="*66)
for metric in ["ESRRB","mTOR","sterol"]:
    c=results["Ctrl"][metric][results["Ctrl"]["ct"]=="EPI"]
    r=results["RapaLink"][metric][results["RapaLink"]["ct"]=="EPI"]
    c=c[~np.isnan(c)]; r=r[~np.isnan(r)]
    _,p=mannwhitneyu(r,c)
    hyp="(pre-naive model predicts RapaLink HIGHER for ESRRB, LOWER for mTOR/sterol)" if metric!="ESRRB" else "(model predicts RapaLink HIGHER)"
    print(f"  {metric:7s}: Ctrl median={np.median(c):+.2f}, RapaLink median={np.median(r):+.2f}, p={p:.2e}  {hyp}")
print("\nDONE")
