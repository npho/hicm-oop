# ============================================================================
# Figure S5  --  figure build
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   The sterol programme tracks mTOR activity in mouse paused ESCs and
#   diapause - links the axis to a known mTOR-low state.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   figure build
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S5. The sterol programme tracks mTOR activity in mouse paused
#   ESCs and diapause (secondary support). Left: sterol delta-score is lower
#   in mTOR-inhibition-paused ESCs than in active (2i/serum) ESCs across
#   three lines. Right: sterol delta-score is lower in diapaused than
#   reactivated mouse embryos. Both are low-mTOR states, consistent with
#   mTOR acting upstream of the sterol programme.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: name)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/mouse_sterol_fig.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np
fig,axes=plt.subplots(1,2,figsize=(12,5))
# paused ESC
lines=["E14","Oct4-GiP","v6.5"]
base=[np.mean([4.79,3.82]),np.mean([4.56,3.57]),np.mean([4.55,4.04])]
paused=[2.67,3.21,2.86]
x=np.arange(3); w=0.35
axes[0].bar(x-w/2,base,w,label="2i/Serum (mTOR active)",color="#3b7dd8")
axes[0].bar(x+w/2,paused,w,label="Paused (mTOR inhibited)",color="#dc2626")
axes[0].set_xticks(x); axes[0].set_xticklabels(lines); axes[0].set_ylabel("sterol program delta-score")
axes[0].set_title("Mouse paused ESCs (GSE81285)\nsterol LOWER in paused (3/3 lines)"); axes[0].legend(fontsize=8)
# diapause
axes[1].bar([0,1],[2.78,4.03],color=["#dc2626","#3b7dd8"],width=0.6)
axes[1].set_xticks([0,1]); axes[1].set_xticklabels(["Diapause\n(mTOR-low)","Reactivated\n(mTOR active)"])
axes[1].set_ylabel("sterol program delta-score")
axes[1].set_title("Mouse diapause (GSE241462)\nsterol LOWER in diapause")
plt.suptitle("The sterol programme tracks mTOR activity across species and states",fontweight="bold",fontsize=13)
plt.tight_layout(); plt.savefig("/Users/adelosan/cesfw_runs/manuscript_build/figures/fig17_mouse_sterol.png",dpi=150,bbox_inches="tight")
print("saved fig17_mouse_sterol.png")
