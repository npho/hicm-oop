# ============================================================================
# Figure S17  --  TOBIAS macOS fork-mode wrapper
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   ICM versus naive-hESC (epiblast proxy) accessibility.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   TOBIAS macOS fork-mode wrapper
#
# PUBLISHED CAPTION (verbatim, for cross-checking):
#   Figure S17. ICM-versus-naïve-hESC (epiblast proxy) accessibility. The
#   confident, promoter-anchored normal-ploidy ICM peaks are shared with the
#   naïve/epiblast-like state, consistent with a shared (not ICM-specific)
#   circuit. No human embryo epiblast ATAC exists; naïve hESC is used as the
#   epiblast proxy.
#
# ATTRIBUTION CONFIDENCE: HIGH   (evidence: ms)
#   hash = embedded manuscript image is byte-identical to a PNG this
#          script writes; ms = named in the manuscript Methods;
#          name = output filename matches unambiguously;
#          infer = attributed from caption + script contents, VERIFY.
# ORIGINAL PATH: /Users/adelosan/cESFW_analysis/atac_rt_pipeline/tobias_fork.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

#!/usr/bin/env python
"""Force multiprocessing 'fork' start method so TOBIAS's mp.Pool does not try to
pickle its logger (macOS Python 3.8+ defaults to 'spawn', which crashes TOBIAS
with PicklingError: logger cannot be pickled). Then dispatch to the TOBIAS CLI."""
import multiprocessing as mp
mp.set_start_method("fork", force=True)
from tobias.TOBIAS import main
main()
