# ============================================================================
# Figure S1  --  underlying statistics
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Branch/ICM cells are not stressed or apoptotic 'reject' cells - rules
#   out the main alternative explanation for a transcriptionally
#   intermediate population.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   underlying statistics
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S1. cESFW branch/ICM cells show no distinct stress/apoptosis
#   “reject” signature. Apoptosis and stress (unfolded-protein, heat-shock
#   and oxidative) programme delta-scores, percentage mitochondrial reads,
#   and genes detected per cell across human embryo lineages in
#   developmental order; the questioned cESFW groups (ICM/TE Branch, ICM,
#   Epi/Hypo Branch; red) versus other lineages (grey). Apoptosis scores in
#   the questioned groups lie within the range of bona fide stages (8-cell,
#   morula and embryonic disc are comparable or higher); stress programmes
#   are not elevated (if anything lower); the modestly lower complexity of
#   the ICM and Epi/Hypo Branch is consistent with the low transcriptional
#   output of the ICM and largely reflects library size. The cells therefore
#   lack the coordinated high-stress, high-apoptosis, low-quality profile
#
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_reject_cells.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Are the cESFW branch/ICM cells "REJECT" cells (stressed / apoptotic / low-quality),
or genuine transitional states?  Tests stress & apoptosis programs + QC metrics
per lineage in the human embryo reference, comparing the questioned groups
(ICM/TE Branch, ICM, Epi/Hypo Branch) with well-defined preimplantation lineages.

Run: conda run -n cesfw_env python analysis_reject_cells.py
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _load_data import load_embryo_data
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")

APOPTOSIS = ["BAX","BAK1","BID","BBC3","PMAIP1","BCL2L11","CASP3","CASP6","CASP7","CASP8",
             "CASP9","APAF1","CYCS","FAS","FADD","TNFRSF10B","DIABLO","TP53","CDKN1A","GADD45A"]
STRESS = ["DDIT3","ATF4","ATF3","XBP1","HSPA5","HERPUD1","DNAJB9","PPP1R15A",   # UPR
          "HSPA1B","HSPA6","HSP90AA1","HSP90AB1","DNAJB1","DNAJA1","HSPB1","BAG3",  # heat-shock
          "HMOX1","NQO1","SOD2","GPX1","GCLM","TXN","TXNRD1"]                   # oxidative

SI, C, genes, E, *_ = load_embryo_data()
ann = SI["Manual_Annotations"].astype(str).values
mt = [g for g in C.columns if g.upper().startswith("MT-")]

ad = sc.AnnData(X=C.values.astype(float)); ad.var_names = C.columns; ad.obs_names = C.index
ad.obs["ann"] = ann
ad.obs["total_counts"] = C.values.sum(1)
ad.obs["n_genes"] = (C.values > 0).sum(1)
ad.obs["pct_mito"] = 100 * C[mt].values.sum(1) / np.maximum(ad.obs["total_counts"].values, 1)
sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
sc.tl.score_genes(ad, [g for g in APOPTOSIS if g in ad.var_names], score_name="apoptosis", ctrl_size=50)
sc.tl.score_genes(ad, [g for g in STRESS if g in ad.var_names], score_name="stress", ctrl_size=50)

QUESTIONED = ["ICM/TE Branch", "ICM", "Epi/Hypo Branch"]
# well-defined preimplantation reference lineages (not branch, not Unknown)
REFERENCE = ["8-Cell", "Morula", "preIm-Epi", "Hypo", "Early TE", "Mid TE", "Polar TE", "Mural TE"]

df = ad.obs.copy()
print("=== Per-lineage medians ===")
tab = df.groupby("ann")[["apoptosis","stress","pct_mito","n_genes","total_counts"]].median()
order = ["8-Cell","Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo",
         "Embryonic Disc","Early TE","Mid TE","Polar TE","Mural TE","cTB","sTB","ExE-Mech","Unknown"]
tab = tab.reindex([o for o in order if o in tab.index])
print(tab.round(3).to_string())
tab.to_csv(os.path.join(OUT, "reject_per_lineage_medians.csv"))

ref_mask = df["ann"].isin(REFERENCE)
print("\n=== Questioned group vs pooled well-defined preimplantation reference ===")
print(f"(reference = {REFERENCE}, n={ref_mask.sum()})")
rows=[]
for grp in QUESTIONED:
    gm = df["ann"] == grp
    for metric, worse in [("apoptosis","higher"),("stress","higher"),
                          ("pct_mito","higher"),("n_genes","lower")]:
        a = df.loc[gm, metric].values; b = df.loc[ref_mask, metric].values
        U,p = mannwhitneyu(a,b,alternative="two-sided"); rbc = 1-2*U/(len(a)*len(b))
        rows.append(dict(group=grp, metric=metric, grp_median=np.median(a),
                         ref_median=np.median(b), p=p, effect_rbc=rbc))
        flag = ""
        if metric in ("apoptosis","stress","pct_mito") and np.median(a)>np.median(b) and p<0.05: flag=" <-- elevated"
        if metric=="n_genes" and np.median(a)<np.median(b) and p<0.05: flag=" <-- reduced"
        print(f"  {grp:16s} {metric:12s}: grp={np.median(a):8.3f} ref={np.median(b):8.3f} "
              f"| p={p:.1e} rbc={rbc:+.2f}{flag}")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "reject_questioned_vs_reference.csv"), index=False)
df[["ann","apoptosis","stress","pct_mito","n_genes","total_counts"]].to_csv(os.path.join(OUT,"reject_cell_scores.csv"))

# ── library-size / complexity ADJUSTED test (does the signal survive depth?) ──
import statsmodels.formula.api as smf
df["log10_counts"] = np.log10(df["total_counts"] + 1)
sub = df[df["ann"].isin(QUESTIONED + REFERENCE)].copy()
print("\n=== Library-size + complexity ADJUSTED (OLS: score ~ group + log10_counts + n_genes) ===")
adj=[]
for grp in QUESTIONED:
    d = sub[sub["ann"].isin([grp] + REFERENCE)].copy(); d["q"] = (d["ann"] == grp).astype(int)
    for metric in ["apoptosis","stress"]:
        m = smf.ols(f"{metric} ~ q + log10_counts + n_genes", data=d).fit()
        b, p = m.params["q"], m.pvalues["q"]
        adj.append(dict(group=grp, metric=metric, adj_beta=b, adj_p=p))
        print(f"  {grp:16s} {metric:10s}: adj beta={b:+.3f} p={p:.1e} "
              f"{'(elevated)' if b>0 and p<0.05 else '(not elevated)'}")
# complexity relative to library size (depth-adjusted residual gene count)
base = smf.ols("n_genes ~ log10_counts", data=sub).fit()
sub["resid_ngenes"] = sub["n_genes"] - base.predict(sub)
print("  depth-adjusted residual genes/cell:",
      {g: int(sub[sub['ann']==g]['resid_ngenes'].median()) for g in QUESTIONED})
pd.DataFrame(adj).to_csv(os.path.join(OUT, "reject_adjusted_ols.csv"), index=False)
print("\nSaved -> reject_* in", OUT)
