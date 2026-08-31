# ============================================================================
# Figure S19  --  per-study replication of programme shapes
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Inter-study robustness of every key ICM analysis, ordered by main
#   figure.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   per-study replication of programme shapes
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
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/transition_replicate.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
from scipy.stats import mannwhitneyu, kruskal

OUT = "/Users/adelosan/cesfw_runs/transition/"
Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
all_genes = list(Counts.columns)

# Load Petropoulos-derived gene programs
petro_shapes = pd.read_csv(OUT + "transition_gene_shapes.csv", index_col=0)
programs = {
    "monotonic_up (biosynthesis/mTORC1/lineage)": petro_shapes[petro_shapes["shape"]=="monotonic_up"].index.tolist(),
    "monotonic_down (ribosome/proliferation)":     petro_shapes[petro_shapes["shape"]=="monotonic_down"].index.tolist(),
    "peak_at_ICM (metallothionein/p53)":           petro_shapes[petro_shapes["shape"]=="peak_at_ICM"].index.tolist(),
}
for k,v in programs.items():
    print(f"{k}: {len(v)} genes")

STAGES = ["ICM/TE Branch", "ICM", "Epi/Hypo Branch"]

def stage_score(dataset, gene_set):
    """Mean z-scored expression of a gene set at each of the 3 stages, in a dataset."""
    d = Sample_Info[Sample_Info.Dataset == dataset]
    cells = {s: d[d.Manual_Annotations == s].index for s in STAGES}
    ns = {s: len(cells[s]) for s in STAGES}
    if min(ns.values()) < 4:
        return None, ns
    allc = np.concatenate([cells[s] for s in STAGES])
    counts_sub = Counts.loc[allc]
    lib = counts_sub.sum(axis=1); lib[lib==0]=1
    log2cpm = np.log2(counts_sub.div(lib, axis=0)*1e6 + 1)
    genes_present = [g for g in gene_set if g in log2cpm.columns]
    # per-cell mean of the program, then per-stage mean
    prog = log2cpm[genes_present].mean(axis=1)
    return {s: prog.loc[cells[s]].mean() for s in STAGES}, ns

print("\n" + "="*70)
print("REPLICATION: do the Petropoulos gene PROGRAMS show the same shape elsewhere?")
print("="*70)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (pname, genes) in zip(axes, programs.items()):
    for ds, color in [("Petropoulos 2016","#16a34a"),("Yanagida 2020","#3b7dd8"),("Meistermann 2021","#f59e0b")]:
        sc, ns = stage_score(ds, genes)
        if sc is None:
            print(f"  {pname[:25]} | {ds}: insufficient cells {ns}")
            continue
        vals = [sc[s] for s in STAGES]
        # normalize to start at 0 for shape comparison
        vals0 = [v - vals[0] for v in vals]
        ax.plot(range(3), vals0, marker="o", label=f"{ds}", color=color, lw=2)
        print(f"  {pname[:35]:35s} | {ds:17s}: {STAGES[0][:6]}={vals[0]:.2f} ICM={vals[1]:.2f} {STAGES[2][:6]}={vals[2]:.2f}  (n={ns})")
    ax.set_xticks(range(3)); ax.set_xticklabels(["ICM/TE","ICM","Epi/Hypo"], fontsize=9)
    ax.set_ylabel("program mean (Δ from ICM/TE)")
    ax.set_title(pname, fontsize=9)
    ax.axhline(0, color="grey", lw=0.5)
    ax.legend(fontsize=7)
plt.suptitle("Do Petropoulos-defined transition programs replicate in Yanagida & Meistermann?", fontweight="bold")
plt.tight_layout()
plt.savefig(OUT + "transition_replication.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: {OUT}transition_replication.png")

# ── Dissect the peak_at_ICM (metallothionein/p53) set specifically ─────────
print("\n" + "="*70)
print("DISSECTING THE ICM-TRANSIENT PEAK (metallothionein/p53)")
print("="*70)
peak_genes = programs["peak_at_ICM (metallothionein/p53)"]
mt_genes = [g for g in peak_genes if g.startswith("MT1") or g.startswith("MT2")]
print(f"Metallothionein genes in the peak set: {mt_genes}")

# Is the MT peak real and specific? Test MT score: ICM vs ICM/TE AND ICM vs Epi/Hypo, per dataset
mt_all = [g for g in all_genes if (g.startswith("MT1") or g.startswith("MT2")) and len(g)<=5]
print(f"\nAll metallothionein genes in data: {mt_all}")
for ds in ["Petropoulos 2016","Yanagida 2020","Meistermann 2021"]:
    d = Sample_Info[Sample_Info.Dataset == ds]
    cells = {s: d[d.Manual_Annotations == s].index for s in STAGES}
    if min(len(cells[s]) for s in STAGES) < 4:
        print(f"  {ds}: insufficient cells"); continue
    allc = np.concatenate([cells[s] for s in STAGES])
    counts_sub = Counts.loc[allc]; lib = counts_sub.sum(axis=1); lib[lib==0]=1
    log2cpm = np.log2(counts_sub.div(lib, axis=0)*1e6 + 1)
    present = [g for g in mt_all if g in log2cpm.columns]
    mt_score = log2cpm[present].mean(axis=1)
    icm_v = mt_score.loc[cells["ICM"]]
    try:
        _, p_lo = mannwhitneyu(icm_v, mt_score.loc[cells["ICM/TE Branch"]], alternative="greater")
        _, p_hi = mannwhitneyu(icm_v, mt_score.loc[cells["Epi/Hypo Branch"]], alternative="greater")
        print(f"  {ds}: MT score ICM/TE={mt_score.loc[cells['ICM/TE Branch']].median():.2f} "
              f"ICM={icm_v.median():.2f} Epi/Hypo={mt_score.loc[cells['Epi/Hypo Branch']].median():.2f} | "
              f"ICM>ICM/TE p={p_lo:.3f}, ICM>Epi/Hypo p={p_hi:.3f}")
    except ValueError as e:
        print(f"  {ds}: test failed {e}")
