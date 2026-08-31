# ============================================================================
# Figure S19  --  per-study TF-activity replication
# ============================================================================
# MANUSCRIPT: "The human inner cell mass is the origin of pluripotency"
#
# WHAT THIS FIGURE SHOWS AND WHY IT IS IN THE PAPER:
#   Inter-study robustness of every key ICM analysis, ordered by main
#   figure.
#
# THIS SCRIPT'S ROLE IN THE FIGURE:
#   per-study TF-activity replication
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
# ORIGINAL PATH: /Users/adelosan/cesfw_runs/tf_replicate.py
# ENVIRONMENT: conda env `cesfw_env`, Python 3.10 (see README.md)
# ============================================================================

import numpy as np, pandas as pd, sys, os
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "/Users/adelosan/cESFW_analysis")
import _load_data_yanagida as ldy
import decoupler as dc
from scipy.stats import spearmanr

Sample_Info, Counts, cESFW_Genes, Embedding, _, _ = ldy.load_embryo_data()
net = pd.read_csv("/Users/adelosan/cesfw_runs/collectri_net.csv")
STAGES = ["Morula","ICM/TE Branch","ICM","Epi/Hypo Branch","preIm-Epi","Hypo"]
focus_tfs = ["SREBF1","SREBF2","PPARA","NFE2L2","TEAD4","SMAD1","SREBF1","AR","SP1"]

print("Replicating TF-activity drivers across datasets (Spearman rho vs DPT pseudotime, embryo-level too)\n")
results = {}
for ds, dptfile in [("Petropoulos 2016","petropoulos_with_dpt.csv"),
                    ("Yanagida 2020","yanagida_with_dpt.csv"),
                    ("Meistermann 2021","meistermann_with_dpt.csv")]:
    try:
        dpt = pd.read_csv(f"/Users/adelosan/cesfw_runs/proper_trajectory/{dptfile}", index_col=0)
    except FileNotFoundError:
        print(f"{ds}: no dpt file"); continue
    cs = Counts.loc[dpt.index]; lib=cs.sum(axis=1); lib[lib==0]=1
    lc = np.log2(cs.div(lib,axis=0)*1e6+1)
    acts = dc.mt.ulm(data=lc, net=net, tmin=5)
    tf = acts[0] if isinstance(acts, tuple) else acts
    dpt2 = dpt.loc[tf.index]; t = dpt2["dpt_pseudotime"].values
    print(f"### {ds} (n={len(tf)} cells, {dpt2['Embryo'].nunique()} embryos) ###")
    for gene in ["SREBF1","SREBF2","PPARA"]:
        if gene not in tf.columns:
            print(f"  {gene}: not in regulon set"); continue
        rho,p = spearmanr(t, tf[gene].values)
        # embryo-level
        emb = pd.DataFrame({"emb":dpt2["Embryo"].values,"dpt":t,"act":tf[gene].values}).groupby("emb").mean()
        rho_e,p_e = spearmanr(emb["dpt"], emb["act"]) if len(emb)>=5 else (np.nan,np.nan)
        print(f"  {gene}: cell rho={rho:+.3f} (p={p:.1e}) | embryo rho={rho_e:+.3f} (p={p_e:.3f}, n={len(emb)})")
    print()
