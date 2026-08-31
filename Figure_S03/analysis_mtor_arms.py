# ============================================================================
# Figure S3  --  underlying arm decomposition
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Metabolic specificity of the ICM transition (which programmes switch,
#   and that the effect is not generic to all metabolism).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   underlying arm decomposition
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S3. Metabolic specificity of the ICM transition. (A) Human
#   ICM→preEPI switch magnitude (Cohen’s d) for each mTORC1-associated
#   biosynthetic arm; sterol and one-carbon are largest (cell-cycle/depth-
#   adjusted). (B) Cross-species: mouse E3.5→E4.5 cell-cycle/depth-adjusted
#   coefficient for each arm in two independent single-cell datasets. Sterol
#   biosynthesis is reproducibly conserved in both, whereas the one-carbon
#   rise is not (n.s. in Mohammed), consistent with sterol being the
#   conserved biosynthetic arm. The bulk ribosome-biogenesis programme is
#   unchanged across the pre-implantation transition (main text), indicating
#   a selective sterol-focused upregulation rather than uniform anabolic
#   activation.
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_mtor_arms.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
mTORC1 downstream-arm decomposition across the human ICM->epiblast transition.
Shows WHICH mTORC1 outputs rise vs fall, with embryo-clustered regression
adjusted for cell-cycle + depth. Tests whether cholesterol is the main arm or one
of several, and whether the translational arm goes the OTHER way (selectivity).
Run: conda run -n cesfw_env python analysis_mtor_arms.py
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")

LINE = ["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Embryonic Disc","Hypo"]
STAGE_ORDER = ["ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Embryonic Disc"]

ARMS = {
 "Protein synthesis / ribosome": None,  # filled from RPL/RPS below
 "Cholesterol synthesis (SREBP2)": ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1",
     "SQLE","LSS","CYP51A1","TM7SF2","MSMO1","NSDHL","HSD17B7","SC5D","DHCR7","DHCR24","LDLR","INSIG1","SREBF2"],
 "Fatty-acid synthesis (SREBP1)": ["FASN","ACACA","ACACB","ACLY","SCD","ELOVL6","SREBF1","ACSS2","MID1IP1"],
 "Nucleotide synthesis": ["CAD","DHODH","UMPS","PPAT","ATIC","GART","PAICS","IMPDH1","IMPDH2","ADSL"],
 "One-carbon / serine": ["PHGDH","PSAT1","PSPH","SHMT1","SHMT2","MTHFD1","MTHFD2","MTHFD1L"],
 "Autophagy / lysosomal (TFEB)": ["MAP1LC3B","SQSTM1","ATG5","ATG7","ULK1","GABARAPL1","LAMP1",
     "LAMP2","CTSB","CTSD","TFEB","WIPI2","BECN1","ATG12"],
}
S_G=["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","MCM6","CDCA7","UHRF1","GINS2","RRM2","CDC6","CDC45","RAD51","GMNN","SLBP"]
G2M=["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","CKS2","MKI67","CCNB2","AURKB","BUB1","CDC20","CENPA","CENPF"]

SI, C, g, E, *_ = load_embryo_data()
ann = SI["Manual_Annotations"].astype(str).values
ARMS["Protein synthesis / ribosome"] = [x for x in C.columns if x.startswith(("RPL","RPS")) and "P" not in x[-2:]][:80]

mask = np.isin(ann, LINE)
ad = sc.AnnData(X=C.values[mask].astype(float)); ad.var_names=C.columns; ad.obs_names=C.index[mask]
ad.obs["ann"]=ann[mask]; ad.obs["embryo"]=SI["Embryo"].values[mask].astype(str)
ad.obs["log_depth"]=np.log10(C.values[mask].sum(1)+1)
sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
sc.pp.neighbors(ad,n_neighbors=15,use_rep="X",metric="correlation"); sc.tl.diffmap(ad)
ad.uns["iroot"]=int(np.where(ad.obs["ann"].values=="ICM/TE Branch")[0][0]); sc.tl.dpt(ad)
ad.obs["pt"]=ad.obs["dpt_pseudotime"]
sc.tl.score_genes_cell_cycle(ad,s_genes=[x for x in S_G if x in ad.var_names],g2m_genes=[x for x in G2M if x in ad.var_names])
for nm,genes in ARMS.items():
    sc.tl.score_genes(ad,[x for x in genes if x in ad.var_names],score_name=nm,ctrl_size=50)

d=ad.obs.copy(); d["ec"]=d["embryo"].astype("category").cat.codes
print(f"cells={len(d)}, embryos={d.embryo.nunique()}\n")
print("=== mTORC1-arm decomposition: slope vs pseudotime (embryo-clustered, +S+G2M+depth) ===")
rows=[]
for nm in ARMS:
    m=smf.ols(f'Q("{nm}") ~ pt + S_score + G2M_score + log_depth',data=d).fit(
        cov_type="cluster",cov_kwds={"groups":d["ec"].values})
    b,p=m.params["pt"],m.pvalues["pt"]
    rows.append((nm,b,p)); arrow="UP" if b>0 else "DOWN"
    print(f"  {nm:34s}: beta_pt={b:+.3f}  p={p:.1e}  [{arrow}]")
pd.DataFrame(rows,columns=["arm","beta_pt","p"]).to_csv(os.path.join(OUT,"mtor_arms_summary.csv"),index=False)

# ── figure: stage-ordered mean score per arm ─────────────────────────────────
fig,ax=plt.subplots(1,2,figsize=(15,5.8))
colors=plt.cm.tab10(np.linspace(0,1,len(ARMS)))
for (nm,_,_),col in zip(rows,colors):
    means=[d[d.ann==s][nm].mean() for s in STAGE_ORDER]
    # z-normalise for comparability of shape
    mz=(np.array(means)-np.mean(means))/(np.std(means)+1e-9)
    ax[0].plot(range(len(STAGE_ORDER)),mz,"o-",color=col,label=nm,lw=2)
ax[0].set_xticks(range(len(STAGE_ORDER))); ax[0].set_xticklabels(STAGE_ORDER,rotation=30,ha="right",fontsize=8)
ax[0].set_ylabel("program score (z across stages)"); ax[0].axhline(0,c="grey",lw=.5)
ax[0].set_title("mTORC1 downstream arms across the ICM transition"); ax[0].legend(fontsize=7.5,frameon=False)
b=[r[1] for r in rows]; nms=[r[0] for r in rows]
bar_c=["#D64550" if x>0 else "#4C72B0" for x in b]
ax[1].barh(range(len(nms)),b,color=bar_c); ax[1].set_yticks(range(len(nms))); ax[1].set_yticklabels(nms,fontsize=8.5)
ax[1].axvline(0,c="k",lw=.8); ax[1].set_xlabel("pseudotime slope (embryo-clustered, CC+depth adj)")
ax[1].set_title("Direction & magnitude per arm"); ax[1].invert_yaxis()
for sp in ["top","right"]: ax[1].spines[sp].set_visible(False)
fig.tight_layout(); out=os.path.join(OUT,"mTOR_arms_decomposition.png")
fig.savefig(out,dpi=200,bbox_inches="tight"); print("\nSaved:",out)
