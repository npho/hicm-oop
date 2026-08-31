# ============================================================================
# Figure S3  --  A,B: ICM->preEPI switch magnitude and mTOR arm decomposition
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Metabolic specificity of the ICM transition (which programmes switch,
#   and that the effect is not generic to all metabolism).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   A,B: ICM->preEPI switch magnitude and mTOR arm decomposition
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
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: hash-name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/fig_mtor_arms_supp.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Supplementary Fig S3: metabolic specificity of the ICM transition.
(A) Human ICM->preEPI switch magnitude per mTORC1 arm (Cohen's d).
(B) Cross-species conservation (mouse E3.5->E4.5, cell-cycle/depth-adjusted beta,
    both datasets) -> sterol conserves strongly, one-carbon does not.
Values from: human_arms.py (A), analysis_mouse_arms.py (B), analysis_mtor_arms.py.
"""
import os
import numpy as np, matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots", "revision_results")

# (A) human ICM->preEPI Cohen's d (CC/depth-adjusted p all <1e-3)
arms = ["Sterol\n(SREBP2)", "One-carbon", "Fatty-acid\n(SREBP1)", "Ribosome"]
human_d = [1.09, 1.42, 0.82, 0.83]
# (B) mouse E3.5->E4.5 adjusted beta (Mohammed, Nowotschin); * = n.s. in that dataset
moh = [0.205, 0.046, 0.114, 0.224]
now = [0.315, 0.063, 0.192, 0.010]
moh_ns = [False, True, False, False]; now_ns = [False, False, False, True]

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
c = ["#C1362F", "#E8A33A", "#E0736A", "#8CA0B3"]
ax[0].bar(range(len(arms)), human_d, color=c)
ax[0].set_xticks(range(len(arms))); ax[0].set_xticklabels(arms, fontsize=9)
ax[0].set_ylabel("ICM $\\rightarrow$ preEPI switch (Cohen's d)")
ax[0].set_title("A  Human: switch magnitude\n(all arms rise; one-carbon & sterol largest)", fontsize=11)
ax[0].spines[["top", "right"]].set_visible(False)

x = np.arange(len(arms)); w = 0.38
b1 = ax[1].bar(x - w/2, moh, w, color="#4C72B0", label="Mohammed 2017")
b2 = ax[1].bar(x + w/2, now, w, color="#8FB8DE", label="Nowotschin 2019")
for i, ns in enumerate(moh_ns):
    if ns: ax[1].text(x[i]-w/2, moh[i]+0.01, "n.s.", ha="center", fontsize=8, color="grey")
for i, ns in enumerate(now_ns):
    if ns: ax[1].text(x[i]+w/2, now[i]+0.01, "n.s.", ha="center", fontsize=8, color="grey")
ax[1].set_xticks(x); ax[1].set_xticklabels(arms, fontsize=9)
ax[1].set_ylabel("mouse E3.5$\\rightarrow$E4.5 adj. $\\beta$")
ax[1].set_title("B  Cross-species conservation\n(sterol conserves; one-carbon does not)", fontsize=11)
ax[1].legend(frameon=False, fontsize=9); ax[1].axhline(0, c="k", lw=.6)
ax[1].spines[["top", "right"]].set_visible(False)

fig.tight_layout()
out = os.path.join(OUT, "FigS3_metabolic_specificity.png")
fig.savefig(out, dpi=250, bbox_inches="tight"); print("Saved:", out)
