# ============================================================================
# Figure S21  --  genome-wide ignition/transience scan
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM as the NANOG/GATA6 balance point, and which factors meet the
#   ignition criteria once trophectoderm is used as an out-group.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   genome-wide ignition/transience scan
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S21. The ICM as the point where the NANOG and GATA6 programmes
#   come into balance, and which factors meet the ignition criteria. (A)
#   Mean expression (+/- s.e.m.) of NANOG and GATA6 along the pluripotent
#   lineage; NANOG rises and GATA6 falls, and the two converge at the ICM
#   (shaded), crossing immediately afterwards. (B) NANOG minus GATA6 for
#   every blastocyst compartment. The ICM (red) is the only compartment near
#   zero; the ICM/TE Branch (orange) is still GATA6-dominant, so the balance
#   tips within the ICM-family. (C) The two ignition criteria plotted
#   against each other for all 13,041 detected genes (grey), with
#   trophectoderm restriction given in parentheses for the labelled factors.
#   ESRRB occupies the ignited-and-transient quadrant and is 2.6 log2 units
#   higher in the ICM than across trophectoderm; TFCP2L1 has the larger
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: authored-here)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/ignition_search.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
What actually 'ignites' at the human ICM?
Genome-wide search for genes that satisfy the paper's own ignition criteria along the
PLURIPOTENT lineage (up at ICM vs upstream, down toward epiblast = transient), and then
asks which of those are ALSO not shared with trophectoderm -- the property ESRRB and
TFCP2L1 turn out to lack.
"""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

S, C, cesfw_genes, _, _, _ = ldy.load_embryo_data()
df = S[["Dataset", "Manual_Annotations"]].copy()
lib = C.sum(axis=1); lib[lib == 0] = 1
lc = np.log2(C.div(lib, axis=0) * 1e6 + 1).reindex(df.index)

UP   = ["Morula"]                        # immediately upstream of the ICM-family
ICMS = ["ICM", "ICM/TE Branch"]          # the ICM-family (paper's framework)
DOWN = ["Epi/Hypo Branch", "preIm-Epi"] # downstream, epiblast-directed
TE   = ["Early TE", "Mid TE", "Mural TE", "Polar TE"]

def gmean(stages):
    m = df.Manual_Annotations.isin(stages).values
    return lc[m].mean(axis=0), m.sum()

up, n_up = gmean(UP); icm, n_icm = gmean(ICMS)
dn, n_dn = gmean(DOWN); te, n_te = gmean(TE)
print(f"n: upstream={n_up}, ICM={n_icm}, downstream={n_dn}, TE={n_te}")

det = (lc > 0).mean(axis=0)
keep = det[det > 0.10].index                       # stainable / real, not noise
up, icm, dn, te = up[keep], icm[keep], dn[keep], te[keep]

res = pd.DataFrame({
    "ICM": icm.round(2),
    "vs_upstream": (icm - up).round(2),   # ignition: up at the ICM
    "vs_epiblast": (icm - dn).round(2),   # transient: down again after
    "vs_TE": (icm - te).round(2),         # lineage-restricted: not a TE program
})
# ignition strength = the weaker of the two within-lineage contrasts
res["ignition"] = res[["vs_upstream", "vs_epiblast"]].min(axis=1).round(2)
res["in_cESFW"] = [g in set(cesfw_genes) for g in res.index]

print("\n" + "=" * 104)
print("A) IGNITION ALONG THE PLURIPOTENT LINEAGE ONLY (paper's criterion): top 25")
print("=" * 104)
a = res.sort_values("ignition", ascending=False).head(25)
print(a.to_string())
for g in ["ESRRB", "TFCP2L1", "NANOG", "GATA6", "KLF4", "KLF5", "TBX3"]:
    if g in res.index:
        r = res.loc[g]
        rank = int((res.ignition > r.ignition).sum() + 1)
        print(f"   {g:8s} ignition={r.ignition:6.2f}  rank {rank:5d}/{len(res)}   "
              f"vs_up={r.vs_upstream:6.2f} vs_epi={r.vs_epiblast:6.2f} vs_TE={r.vs_TE:6.2f}")

print("\n" + "=" * 104)
print("B) IGNITION *AND* NOT A TROPHECTODERM PROGRAM (ignition>0 and vs_TE>1): top 25")
print("=" * 104)
b = res[(res.ignition > 0) & (res.vs_TE > 1)].sort_values("ignition", ascending=False)
print(f"n genes passing both = {len(b)}")
print(b.head(25).to_string())

print("\n" + "=" * 104)
print("C) Where ESRRB / TFCP2L1 fail the TE filter")
print("=" * 104)
print(res.loc[[g for g in ["ESRRB", "TFCP2L1"] if g in res.index]].to_string())

print("\n" + "=" * 104)
print("D) Do the B-list genes hold up per dataset? (top 10, ICM minus upstream)")
print("=" * 104)
top10 = list(b.head(10).index)
for ds in ["Petropoulos 2016", "Yanagida 2020", "Meistermann 2021"]:
    dm = (df.Dataset == ds).values
    mu = dm & df.Manual_Annotations.isin(UP).values
    mi = dm & df.Manual_Annotations.isin(ICMS).values
    if mi.sum() < 5 or mu.sum() < 5:
        print(f"{ds}: too few cells (ICM={mi.sum()}, up={mu.sum()})"); continue
    d = (lc.loc[mi, top10].mean() - lc.loc[mu, top10].mean()).round(2)
    print(f"{ds} (ICM n={mi.sum()}, upstream n={mu.sum()}):")
    print("   " + d.to_string().replace("\n", "\n   "))
