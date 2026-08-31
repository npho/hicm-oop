# ============================================================================
# Figure 4  --  H: ESRRB-independent ICM Identity Score for the Dattani intermediate
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   No in-vitro system stabilises the ICM programme, but rare genuine ESRRB+
#   cells arise. Two independent measurements (reference mapping; direct
#   signature scoring).
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   H: ESRRB-independent ICM Identity Score for the Dattani intermediate
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
# ATTRIBUTION CONFIDENCE: MEDIUM   (evidence: infer)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/analysis_dattani_icmlike.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

"""Test Dattani et al. 2024 (Guo lab, Cell Stem Cell) 'transitional ICM-like
state' against our ICM criterion: is it ESRRB+GATA6+ (dual-competent)?

Data: Mendeley bg92x53rk4, adata_n2b27_fax_fig5.h5ad (naive PSC -> hypoblast
FGF-driven differentiation timecourse; the ICM-like state appears mid-course).

Step 1 inspects the annotation columns/labels (we don't know their naming yet).
Step 2 scores ESRRB / GATA6 / co-expression per annotated group.
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

H5 = "/Users/adelosan/cESFW_analysis/data/Dattani2024/adata_n2b27_fax_fig5.h5ad"
ad = sc.read_h5ad(H5)
print("shape:", ad.shape)
print("\n=== obs columns ===")
for c in ad.obs.columns:
    vals = ad.obs[c]
    if vals.dtype.name in ("category", "object") or vals.nunique() < 40:
        print(f"  {c}: {list(pd.Series(vals).astype(str).value_counts().head(20).items())}")

print("\n=== is ESRRB / GATA6 in var? ===")
for g in ["ESRRB", "GATA6", "NANOG", "TFCP2L1", "PDGFRA", "SOX17"]:
    print(f"  {g}: {'yes' if g in ad.var_names else 'NO'}")
print("var_names sample:", list(ad.var_names[:5]))
# layers / raw
print("\nlayers:", list(ad.layers.keys()), "| raw:", ad.raw is not None, "| X max:", float(ad.X.max()))
