# ============================================================================
# Figure S21  --  A-D: crossover, balance, ignition/transience scatter, per-cell modality
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The ICM as the NANOG/GATA6 balance point, and which factors meet the
#   ignition criteria once trophectoderm is used as an out-group.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A-D: crossover, balance, ignition/transience scatter, per-cell modality
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
# ORIGINAL PATH: /private/tmp/claude-501/-Users-adelosan-Downloads-threshold-ptd/b69a1676-2239-4389-9e30-373e2af378bb/scratchpad/make_figS21.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Figure S21: the ICM as the NANOG/GATA6 crossover, and which factors meet the
ignition criteria once trophectoderm is used as an out-group."""
import numpy as np, pandas as pd, sys
import warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy

OUT = "/Users/adelosan/cESFW_analysis/plots/Fig_S21_ICM_crossover_ignition.png"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
                     "font.size": 7, "axes.linewidth": 0.6})

S, C, cesfw_genes, _, _, _ = ldy.load_embryo_data()
df = S[["Dataset", "Manual_Annotations"]].copy()
lib = C.sum(axis=1); lib[lib == 0] = 1
lc = np.log2(C.div(lib, axis=0) * 1e6 + 1).reindex(df.index)

ORD = ["Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi"]
LBL = ["Morula", "ICM/TE\nBranch", "ICM", "Epi/Hypo\nBranch", "preIm-\nEpi"]
OUTG = {"Hypo": "Hypo", "Early TE": "Early TE", "Mid TE": "Mid TE"}

fig = plt.figure(figsize=(10.6, 7.0))
gs = fig.add_gridspec(2, 2, wspace=0.28, hspace=0.55, left=0.075, right=0.985, top=0.90, bottom=0.075)

# ---- A: NANOG / GATA6 crossover along the pluripotent lineage ----
axA = fig.add_subplot(gs[0, 0])
x = np.arange(len(ORD))


def stage_mean(g, s):
    m = (df.Manual_Annotations == s).values
    return lc.loc[m, g].mean()


def stage_sem(g, s):
    m = (df.Manual_Annotations == s).values
    v = lc.loc[m, g]
    return v.std() / np.sqrt(max(len(v), 1))


for g, col, mk in [("NANOG", "#8172B2", "D"), ("GATA6", "#55A868", "^")]:
    mu = np.array([stage_mean(g, s) for s in ORD])
    se = np.array([stage_sem(g, s) for s in ORD])
    axA.plot(x, mu, marker=mk, color=col, lw=1.7, ms=5, label=g, zorder=3)
    axA.fill_between(x, mu - se, mu + se, color=col, alpha=0.18, lw=0, zorder=2)
axA.axvspan(1.5, 2.5, color="#FFD54F", alpha=0.22, zorder=0)
axA.set_xticks(x); axA.set_xticklabels(LBL, fontsize=6.4)
axA.set_ylabel("mean log$_2$CPM (± s.e.m.)")
axA.legend(fontsize=6.5, frameon=False, loc="upper right", handlelength=1.4)
axA.set_title("A   NANOG rises and GATA6 falls across the ICM;\n"
              "      the two profiles converge there (shaded)",
              fontsize=7.7, loc="left", fontweight="bold")
for sp in ("top", "right"):
    axA.spines[sp].set_visible(False)

# ---- B: NANOG - GATA6 balance, all compartments ----
axB = fig.add_subplot(gs[0, 1])
ALL = ["Morula", "ICM/TE Branch", "ICM", "Epi/Hypo Branch", "preIm-Epi",
       "Hypo", "Early TE", "Mid TE", "Mural TE", "Polar TE"]
ALLL = ["Morula", "ICM/TE Br", "ICM", "Epi/Hypo Br", "preIm-Epi",
        "Hypo", "Early TE", "Mid TE", "Mural TE", "Polar TE"]
bal = np.array([stage_mean("NANOG", s) - stage_mean("GATA6", s) for s in ALL])
cols = ["#D62728" if s == "ICM" else ("#F0A05A" if s == "ICM/TE Branch" else "#9AA5B1")
        for s in ALL]
yb = np.arange(len(ALL))[::-1]
axB.barh(yb, bal, color=cols, edgecolor="none", height=0.7)
axB.axvline(0, color="k", lw=0.9)
axB.set_yticks(yb); axB.set_yticklabels(ALLL, fontsize=6.4)
axB.set_xlabel("NANOG − GATA6  (log$_2$CPM)")
for i, v in zip(yb, bal):
    axB.text(v + (0.18 if v >= 0 else -0.18), i, f"{v:+.1f}", va="center",
             ha="left" if v >= 0 else "right", fontsize=5.8)
axB.set_xlim(-6.4, 5.0)
axB.set_title("B   The ICM is the only compartment whose mean\n      lies near the balance point",
              fontsize=7.7, loc="left", fontweight="bold")
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)

# ---- C: ignition criteria with TE as out-group ----
axC = fig.add_subplot(gs[1, 0])
UP, ICMF = ["Morula"], ["ICM", "ICM/TE Branch"]
DOWN = ["Epi/Hypo Branch", "preIm-Epi"]
TE = ["Early TE", "Mid TE", "Mural TE", "Polar TE"]


def gm(stages):
    return lc[df.Manual_Annotations.isin(stages).values].mean(axis=0)


up, icmv, dn, te = gm(UP), gm(ICMF), gm(DOWN), gm(TE)
det = (lc > 0).mean(axis=0)
keep = det[det > 0.10].index
vs_up = (icmv - up)[keep]
vs_epi = (icmv - dn)[keep]
vs_te = (icmv - te)[keep]

axC.scatter(vs_up, vs_epi, s=1.2, c="#D8DEE4", alpha=0.5, lw=0, rasterized=True)
axC.axhline(0, color="#888", lw=0.7, ls=":")
axC.axvline(0, color="#888", lw=0.7, ls=":")
MARK = {"ESRRB": "#D62728", "TFCP2L1": "#4C72B0", "KLF4": "#7F7F7F",
        "GATA6": "#55A868", "NANOG": "#8172B2"}
for g, col in MARK.items():
    if g not in vs_up.index:
        continue
    axC.scatter([vs_up[g]], [vs_epi[g]], s=46, c=col, edgecolor="white", lw=0.7, zorder=5)
    axC.annotate(f"{g}\n(TE {vs_te[g]:+.1f})", (vs_up[g], vs_epi[g]),
                 textcoords="offset points", xytext=(7, 5), fontsize=5.9, color=col,
                 fontweight="bold", zorder=6)
axC.set_xlabel("ignition:  ICM − morula (log$_2$CPM)")
axC.set_ylabel("transience:  ICM − epiblast")
axC.set_xlim(-3.2, 5.2); axC.set_ylim(-4.2, 4.2)
axC.text(3.4, 3.6, "ignited &\ntransient", fontsize=5.9, color="#555", ha="center")
axC.text(3.4, -3.6, "ignited &\nmaintained", fontsize=5.9, color="#555", ha="center")
axC.set_title("C   Only ESRRB is both ignited and transient;\n"
              "      TFCP2L1 is ignited but maintained",
              fontsize=7.7, loc="left", fontweight="bold")
for sp in ("top", "right"):
    axC.spines[sp].set_visible(False)


# ---- D: per-cell distributions (is the balance a mixture?) ----
axD = fig.add_subplot(gs[1, 1])
try:
    import diptest; HAVE_DIP = True
except Exception:
    HAVE_DIP = False
bal_all = lc["NANOG"] - lc["GATA6"]
parts, labs, dips = [], [], []
for s_ in ORD:
    m = (df.Manual_Annotations == s_).values
    v = bal_all[m].values; v = v[np.isfinite(v)]
    parts.append(v); labs.append(s_.replace(" Branch", "\nBranch").replace("preIm-Epi", "preIm-\nEpi"))
    dips.append(diptest.diptest(v)[1] if (HAVE_DIP and len(v) >= 10) else np.nan)
vp = axD.violinplot(parts, positions=np.arange(len(parts)), widths=0.82,
                    showmeans=False, showextrema=False)
for i, b in enumerate(vp["bodies"]):
    b.set_facecolor("#D62728" if ORD[i] == "ICM" else "#9AA5B1")
    b.set_alpha(0.55); b.set_edgecolor("none")
for i, v in enumerate(parts):
    axD.scatter(np.random.default_rng(i).normal(i, 0.055, len(v)), v, s=1.6,
                c="#37474F", alpha=0.45, lw=0, zorder=3)
    axD.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="k", lw=1.3, zorder=4)
axD.axhline(0, color="k", lw=0.9, ls="--")
axD.set_xticks(np.arange(len(parts))); axD.set_xticklabels(labs, fontsize=6.3)
axD.set_ylabel("NANOG − GATA6 per cell (log$_2$CPM)")
if HAVE_DIP and np.isfinite(dips[2]):
    axD.annotate(f"ICM: unimodal\n(dip p = {dips[2]:.2f}, n = {len(parts[2])})",
                 (2, axD.get_ylim()[1]), textcoords="offset points", xytext=(4, -22),
                 fontsize=6.0, color="#D62728", fontweight="bold")
axD.set_title("D   Single cells straddle the balance point in one\n"
              "      broad distribution, not two committed groups",
              fontsize=7.7, loc="left", fontweight="bold")
for sp in ("top", "right"):
    axD.spines[sp].set_visible(False)

fig.savefig(OUT, dpi=330, bbox_inches="tight", facecolor="white")
print("wrote", OUT)

rows = []
for s in ALL:
    m = (df.Manual_Annotations == s).values
    rows.append({"compartment": s, "n": int(m.sum()),
                 "NANOG": round(stage_mean("NANOG", s), 2),
                 "GATA6": round(stage_mean("GATA6", s), 2),
                 "NANOG_minus_GATA6": round(stage_mean("NANOG", s) - stage_mean("GATA6", s), 2)})
pd.DataFrame(rows).to_csv(
    "/Users/adelosan/cESFW_analysis/plots/ICM_NANOG_GATA6_balance.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("\nignition/transience/TE for key factors:")
for g in MARK:
    if g in vs_up.index:
        print(f"  {g:8s} vs_morula={vs_up[g]:+.2f}  vs_epiblast={vs_epi[g]:+.2f}  vs_TE={vs_te[g]:+.2f}")
