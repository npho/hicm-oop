# ============================================================================
# Figure 4  --  A,B: cESFW reference-mapping of naive/reset PSC lines
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   No in-vitro system stabilises the ICM programme, but rare genuine ESRRB+
#   cells arise. Two independent measurements (reference mapping; direct
#   signature scoring).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A,B: cESFW reference-mapping of naive/reset PSC lines
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms+name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/naive_projection.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from sklearn.neighbors import KNeighborsClassifier

OUT="/Users/adelosan/cesfw_runs/naive_proj/"; os.makedirs(OUT,exist_ok=True)
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
Counts=Counts.fillna(0)
feat=[g for g in cESFW_Genes if g in Counts.columns]
STAGES=["Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
allg=list(Counts.columns); D="/Users/adelosan/cESFW_analysis/data"
def l2cpm(cs):
    lib=cs.sum(axis=1); lib[lib==0]=1; return np.log2(cs.div(lib,axis=0)*1e6+1)
ref_cells=Sample_Info.index[Sample_Info["Manual_Annotations"].isin(STAGES)]
ref_full=l2cpm(Counts.loc[ref_cells]); ref_lc=ref_full[feat].fillna(0).values
ref_lab=Sample_Info.loc[ref_cells,"Manual_Annotations"].values
KNN=KNeighborsClassifier(n_neighbors=20,metric="correlation",weights="distance").fit(ref_lc,ref_lab)
rng=np.random.default_rng(0)

def load(name):
    if name=="t2iLGo (Messmer)":
        C=pd.read_csv(f"{D}/Messmet2019_sc_Counts.csv",index_col=0)
        info=pd.read_csv(f"{D}/E-MTAB-6819_Sample_Info.csv",index_col=0)
        i2=info[["Assay Name","Factor Value[phenotype]","Characteristics[post-analysis well quality]"]].drop_duplicates(); i2.index=i2["Assay Name"]
        nv=i2[(i2["Factor Value[phenotype]"]=="naive")&(i2["Characteristics[post-analysis well quality]"]=="pass")].index
        return C.loc[[c for c in nv if c in C.index]]
    if name=="t2iLGo (Polo)": return pd.read_csv(f"{D}/polo2020/Polo2020_Naive_QC_counts.csv",index_col=0)
    if name=="e4CL (Mazid)":  return pd.read_csv(f"{D}/Mazid2022_sc_Counts.csv",index_col=0).T
    if name=="PXGL (Guo)":    return pd.read_csv(f"{D}/Guo2021_sc_Counts.csv",index_col=0)
    if name=="5iLAF (Bi)":    return pd.read_csv(f"{D}/bi2022/Bi2022_5iLAF_nESC_QC_counts.csv",index_col=0)
    if name=="5iLAF (Moya)":  return pd.read_csv(f"{D}/moya2022/Moya2022_QC_counts.csv",index_col=0)
    if name=="hEPSC":         return pd.read_csv(f"{D}/HEPSC_Counts.csv",index_col=0)

gmt={}
for line in open("/Users/adelosan/.cache/gseapy/Enrichr.MSigDB_Hallmark_2020.gmt"):
    pp=line.strip().split("\t"); gmt[pp[0]]=[g for g in pp[2:] if g]
srebf=["HMGCR","HMGCS1","SQLE","LDLR","FDFT1","DHCR7","MVD","MVK","FASN","ACACA","SCD"]
def delta(lc,genes):
    p=[g for g in genes if g in lc.columns]; bg=[g for g in allg if g not in p and g in lc.columns]
    return lc[p].mean(axis=1)-lc[bg].mean(axis=1)

def prep(cdf):
    if len(set(cdf.index)&set(feat))>len(set(cdf.columns)&set(feat)): cdf=cdf.T
    if cdf.shape[0]>800: cdf=cdf.iloc[rng.choice(cdf.shape[0],800,replace=False)]
    lib=cdf.sum(axis=1); lib[lib==0]=1; lc=np.log2(cdf.div(lib,axis=0)*1e6+1)
    X=pd.DataFrame(0.0,index=lc.index,columns=feat)
    pres=[g for g in feat if g in lc.columns]; X[pres]=lc[pres].fillna(0).values
    return X,lc

names=["t2iLGo (Messmer)","t2iLGo (Polo)","e4CL (Mazid)","PXGL (Guo)","5iLAF (Bi)","5iLAF (Moya)","hEPSC"]
# real ICM/early reference metabolic values
re=Sample_Info.loc[ref_cells,"Manual_Annotations"].isin(["ICM/TE Branch","ICM"]).values
ref_vals={s:delta(ref_full,g)[re].median() for s,g in [("mTORC1",gmt["mTORC1 Signaling"]),("sterol",srebf),("ESRRB",["ESRRB"])]}
print("REAL ICM/early: "+" | ".join(f"{k}={v:+.2f}" for k,v in ref_vals.items()))
print("="*80)
rows=[]
for nm in names:
    try: cdf=load(nm)
    except Exception as e: print(f"{nm}: load fail {e}"); continue
    X,lc=prep(cdf)
    pred=KNN.predict(X.values); conf=KNN.predict_proba(X.values).max(axis=1)
    vc=pd.Series(pred).value_counts(normalize=True)
    icm=np.isin(pred,["ICM/TE Branch","ICM"])
    dist=" ".join(f"{s.split('/')[0][:5]}:{100*vc.get(s,0):.0f}%" for s in STAGES if vc.get(s,0)>=0.05)
    met=""
    if icm.sum()>=3:
        met=" | ICM-mapped metab: "+", ".join(f"{k}={delta(lc,g)[icm].median():+.2f}" for k,g in [("mTOR",gmt["mTORC1 Signaling"]),("sterol",srebf),("ESRRB",["ESRRB"])])
    print(f"\n{nm} (n={X.shape[0]}, conf={conf.mean():.2f})")
    print(f"  stages: {dist}")
    print(f"  ICM/early mapped: {100*icm.mean():.0f}%{met}")
    rows.append({"condition":nm,"n":X.shape[0],"top":vc.idxmax(),"pct_top":round(100*vc.max(),1),
                 "pct_ICMearly":round(100*icm.mean(),1),"conf":round(conf.mean(),2)})
pd.DataFrame(rows).to_csv(OUT+"projection_summary.csv",index=False)
print("\n"+"="*80); print(pd.DataFrame(rows).to_string(index=False))
print("\nDONE")
