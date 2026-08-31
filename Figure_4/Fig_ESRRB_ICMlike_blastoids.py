# ============================================================================
# Figure 4  --  C,D,E: ICM% + ESRRB+GATA6+ signature + per-cell ESRRB level
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   No in-vitro system stabilises the ICM programme, but rare genuine ESRRB+
#   cells arise. Two independent measurements (reference mapping; direct
#   signature scoring).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   C,D,E: ICM% + ESRRB+GATA6+ signature + per-cell ESRRB level
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
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/Fig_ESRRB_ICMlike_blastoids.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""
Revised Figure 4C: no in-vitro system builds the ICM-like state, shown two ways.

Two independent measurements, plus the ESRRB level in the rare positive cells.

A  ICM by reference-mapping (the cESFW classifier, manuscript Figure 4A/B):
   % of cells assigned ICM identity. Embryo positive control 22.8 %; every
   in-vitro system (naive PSC, blastoids) <= 2.5 %.

B  ICM signature scored directly: % of cells co-expressing the dual-competent
   markers ESRRB+GATA6+ (core pair) and ESRRB+TFCP2L1+GATA6+ (triple). Positive
   control = cESFW embryo ICM-family (ICM + ICM/TE Branch), 93 %; every in-vitro
   system < 5 %. The two measurements agree.

C  How much ESRRB the rare positive cells carry (one dot per ESRRB+ cell, raw
   molecules, log). Kagawa/Yanagida contain a few genuine embryo-level cells
   (156 / 271 molecules, above the 13-molecule ESRRB-negative-control ceiling);
   Iyer never exceeds 2 molecules in Ctrl or RapaLink (mTOR-i); naive/reset PSCs
   and Silva de-diff shown for completeness.

Inputs : plots/ESRRB_ICMsig_dualcompetent.csv, plots/ESRRB_blastoid_allpositives.csv,
         plots/ESRRB_embryo_ICM_positives.csv
Output : plots/Fig_ESRRB_ICMlike_blastoids.{png,pdf}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PLOTS = "/Users/adelosan/cESFW_analysis/plots"
BLUE, ORANGE, VIOLET = "#2a78d6", "#eb6834", "#4a3aa7"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8880", "#e4e3de"
CEIL = 13

sig = pd.read_csv(os.path.join(PLOTS, "ESRRB_ICMsig_dualcompetent.csv"))
pos = pd.read_csv(os.path.join(PLOTS, "ESRRB_blastoid_allpositives.csv"))
emb = pd.read_csv(os.path.join(PLOTS, "ESRRB_embryo_ICM_positives.csv"))
sig = sig.rename(columns={"ICM sig: ESRRB+TFCP2L1+GATA6+%": "triple"})

# systems shared by panels A and B (aligned rows). Classifier % from the
# manuscript reference-mapping; signature % and triple from direct scoring.
# (grp_key for signature lookup, label, classifier%, classifier_txt, double%, is_control)
SYS = [
    ("ICM-family (cESFW)", "Embryo ICM\n(positive control)", 22.8, "22.8%", 93.4, True),
    ("naive hPSC",         "naive hPSC",                     0.5,  "≤0.5%", 0.0,  False),
    ("Yanagida",           "Yanagida 2021",                  0.0,  "0%",    4.5,  False),
    ("Kagawa",             "Kagawa 2022",                    1.0,  "~1%",   0.7,  False),
    ("Iyer (Ctrl)",        "Iyer 2024 (Ctrl)",               2.5,  "2.5%",  2.52, False),
    ("Iyer (RapaLink)",    "Iyer 2024 (RapaLink)",           None, "n.d.",  0.59, False),
]
TRIPLE = {"ICM-family (cESFW)": 92.9, "naive hPSC": 0.0, "Yanagida": 4.5,
          "Kagawa": 0.6, "Iyer (Ctrl)": 2.0, "Iyer (RapaLink)": 0.5}

fig = plt.figure(figsize=(13.6, 5.2), dpi=300)
fig.patch.set_facecolor(SURFACE)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], wspace=0.42,
                      left=0.108, right=0.99, top=0.88, bottom=0.14)
axA, axB, axC = [fig.add_subplot(gs[0, i]) for i in range(3)]
y = np.arange(len(SYS))[::-1]


def col(is_ctrl):
    return BLUE if is_ctrl else ORANGE


# ── A: classifier ────────────────────────────────────────────────────────────
axA.set_facecolor(SURFACE)
for (k, lab, cls, ctxt, dbl, isc), yy in zip(SYS, y):
    if cls is None:
        axA.text(1, yy, "not scored separately", va="center", fontsize=6.6,
                 color=INK_MUTED, style="italic")
        continue
    axA.barh(yy, cls, height=0.6, color=col(isc), zorder=3)
    axA.text(cls + 1.2, yy, ctxt, va="center", fontsize=7.6, color=INK,
             fontweight="semibold")
axA.set_yticks(y)
axA.set_yticklabels([s[1] for s in SYS], fontsize=8, color=INK_2)
axA.set_xlim(0, 30)
axA.set_xticks([0, 10, 20, 30])
axA.set_xticklabels(["0%", "10%", "20%", "30%"], fontsize=8, color=INK_2)
axA.set_xlabel("cells assigned ICM identity", fontsize=8.5, color=INK_2, labelpad=6)
axA.set_title("C   Measurement 1: ICM by\n     cESFW reference-mapping", loc="left",
              fontsize=10, fontweight="semibold", color=INK, pad=8)

# ── B: direct signature ──────────────────────────────────────────────────────
axB.set_facecolor(SURFACE)
h = 0.36
for (k, lab, cls, ctxt, dbl, isc), yy in zip(SYS, y):
    tri = TRIPLE[k]
    axB.barh(yy + h/2, dbl, height=h, color=col(isc), zorder=3)
    axB.barh(yy - h/2, tri, height=h, color=col(isc), alpha=0.5, zorder=3,
             hatch="////", edgecolor=SURFACE, linewidth=0)
    axB.text(dbl + 1.4, yy + h/2, f"{dbl:.1f}%", va="center", fontsize=7.2,
             color=INK, fontweight="semibold")
    axB.text(tri + 1.4, yy - h/2, f"{tri:.1f}%", va="center", fontsize=7.0, color=INK_2)
axB.set_yticks(y)
axB.set_yticklabels([])
axB.set_xlim(0, 112)
axB.set_xticks([0, 25, 50, 75, 100])
axB.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=8, color=INK_2)
axB.set_xlabel("cells carrying the ICM signature", fontsize=8.5, color=INK_2, labelpad=6)
axB.set_title("D   Measurement 2: direct\n     ESRRB⁺GATA6⁺ signature", loc="left",
              fontsize=10, fontweight="semibold", color=INK, pad=8)
axB.legend(handles=[
    Patch(facecolor=INK_MUTED, label="ESRRB⁺ GATA6⁺  (core pair)"),
    Patch(facecolor=INK_MUTED, alpha=0.5, hatch="////", label="+ TFCP2L1⁺  (triple)")],
    loc="lower right", frameon=False, fontsize=6.8, labelcolor=INK_2,
    handlelength=1.3, borderpad=0.2)

# ── C: ESRRB level in positive cells ─────────────────────────────────────────
axC.set_facecolor(SURFACE)
axC.set_xscale("log")
Crows = [("Embryo ICM-family", emb, BLUE),
         ("Yanagida", pos[pos.dataset == "Yanagida"], ORANGE),
         ("Kagawa", pos[pos.dataset == "Kagawa"], ORANGE),
         ("Iyer (Ctrl)", pos[pos.dataset == "Iyer (Ctrl)"], ORANGE),
         ("Iyer (RapaLink)", pos[pos.dataset == "Iyer (RapaLink)"], ORANGE),
         ("naive hPSC (H9)", pos[pos.dataset == "naive hPSC (H9)"], ORANGE),
         ("Silva (de-diff)", pos[pos.dataset == "Silva"], ORANGE)]
yC = np.arange(len(Crows))[::-1]
rng = np.random.default_rng(0)
axC.axvspan(0.6, CEIL, color=VIOLET, alpha=0.06, zorder=0)
axC.axvline(CEIL, color=VIOLET, lw=1.3, ls=(0, (4, 3)), zorder=2)
for (name, sub, c), yy in zip(Crows, yC):
    m = sub.molecules.values
    if len(m) == 0:
        continue
    axC.scatter(m, yy + rng.uniform(-0.22, 0.22, len(m)), s=15, color=c, alpha=0.55,
                edgecolors=SURFACE, linewidths=0.3, zorder=3)
    med = np.median(m)
    axC.plot([med, med], [yy - 0.32, yy + 0.32], color=INK, lw=1.7, zorder=5)
    tag = f"n={len(m)} · max {int(m.max())}"
    if int(m.max()) <= CEIL:
        tag += " ⚠"
    axC.text(6800, yy, tag, va="center", ha="right", fontsize=6.6,
             color=ORANGE if int(m.max()) <= CEIL else INK_MUTED)
axC.text(CEIL * 0.8, len(Crows) - 0.35, f"contamination ceiling\n{CEIL} mol (neg. control)",
         fontsize=6.6, color=VIOLET, ha="right", va="top", linespacing=1.4)
axC.set_yticks(yC)
axC.set_yticklabels([c[0] for c in Crows], fontsize=8, color=INK_2)
axC.set_xlim(0.6, 8000)
axC.set_xticks([1, 10, 100, 1000])
axC.set_xticklabels(["1", "10", "100", "1,000"], fontsize=8, color=INK_2)
axC.set_xlabel("ESRRB molecules per positive cell (log)  ·  ▮ median",
               fontsize=8.5, color=INK_2, labelpad=6)
axC.set_title("E   ESRRB level in the rare positive cells", loc="left",
              fontsize=10, fontweight="semibold", color=INK, pad=8)

for ax in (axA, axB, axC):
    ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, pad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)

# (no suptitle: this renders as the lower half of a merged Figure 4; the
# narrative lives in the figure caption.)

for ext in ("png", "pdf"):
    p = os.path.join(PLOTS, f"Fig_ESRRB_ICMlike_blastoids.{ext}")
    fig.savefig(p, facecolor=SURFACE)
    print("wrote", p)
