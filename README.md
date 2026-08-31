# Figure script audit package
### "The human inner cell mass is the origin of pluripotency"

Generated 2026-08-28 from manuscript `icm_experiment_v57.docx`.

This package contains, for every main and supplementary figure, the script(s) that
produced it. One folder per figure. Every script carries a header block stating what
the figure shows, why it is in the paper, the script's role, the verbatim published
caption, and **how confident the attribution is**.

---

## How each script was attributed to a figure

Attribution was not taken on trust from any existing README (the previous
`paper_code_audit/README.md` maps to an **older figure numbering** and an older title,
and is superseded by this package). Evidence codes appear in every script header:

| code | meaning | strength |
|---|---|---|
| `hash` | the image embedded in the .docx is **byte-identical** (md5) to a PNG on disk, and that PNG filename appears in the script | strongest |
| `ms` | the manuscript Methods **names the script** for this analysis | strong |
| `name` | the script's output filename matches the figure subject unambiguously | good |
| `infer` | attributed from caption content plus reading the script; **auditor should verify** | weakest |

Ground truth for `hash` came from md5-matching all 34 images embedded in the manuscript
against 1,374 candidate PNGs on disk.

---

## Environment

- Python 3.10.14, conda env `cesfw_env`
- scanpy 1.9.8, anndata 0.8.0, scikit-learn 1.2.2, scipy 1.10.1, numpy 1.23.5,
  pandas 1.5.3, statsmodels 0.14.5, decoupler 2.1.4, umap-learn 0.5.5, gseapy 1.1.9,
  harmonypy 2.0.0, diptest 0.11.0, leidenalg 0.12.0, python-igraph 1.0.0,
  matplotlib 3.7.3, cESFW 0.0.1
- Scripts expect data under `/Users/adelosan/cESFW_analysis/data/` and import the
  loaders in `_shared_data_loaders/` — adjust `sys.path` / data paths to relocate.

## Primary data (all public)

| Dataset | Accession |
|---|---|
| Human embryo atlas (cESFW) | Radley et al. 2023, Mendeley 34td4ds2r9 |
| Petropoulos 2016 | E-MTAB-3929 |
| Yanagida 2020/2021 | GSE171820 |
| Meistermann 2021 | in cESFW atlas |
| Marmoset (Boroviak 2018) | E-MTAB-7078 |
| Mouse (Nowotschin 2019) | GSE123046 |
| Mouse diapause | GSE241462 |
| Mouse paused ESC | GSE81285 |
| Iyer 2024 blastoid | GSE267302 |
| Human ICM ATAC (Wu 2018) | GSE101571 |
| Mouse ESC ChIP (Chen 2008) | GSE11431 |
| Braun 2023 brain / Burclaff 2022 gut / Arutyunyan 2023 placenta / Sikkema 2023 + Habermann 2020 lung / Zeng 2025 marrow | see Methods |

---

## Data — where the analysed files came from

**`DATA_PROVENANCE.md`** is the authoritative record: every dataset, its accession, which
files, which figures use it, and where to download it. **`fetch_data.sh`** scripts the
GEO/ArrayExpress downloads (`bash fetch_data.sh /path/to/data`); portal-only datasets
(Mendeley, CELLxGENE, UCSC Cell Browser) are listed for manual retrieval.
**`data_inputs_by_figure.txt`** lists every input file each script reads.

> **Important:** the scripts contain absolute paths of the form
> `/Users/adelosan/cESFW_analysis/data/...` and will not run elsewhere unchanged. Most read
> the embryo reference through `_shared_data_loaders/_load_data_yanagida.py`, whose `DATA`
> constant resolves to `<repo>/data` — repoint that one path and the majority work.

## Composite figures — see `RECONSTRUCTIONS/`

Six figures are composites whose *compositing* step was not found on disk (the component
panels are all present and attributed). `RECONSTRUCTIONS/` contains clearly-labelled
reconstructions that are **not original code**, together with a verification table showing
how closely each matches the published panel.

## Contents

### Figure 1  (`Figure_1/`)

| script | role | confidence |
|---|---|---|
| `radley_original.py` | A: cells on the published cESFW embedding | high (ms+name) |
| `proper_trajectory.py` | B,C: PAGA graph + diffusion pseudotime; writes the *_with_dpt.csv used downstream | high (ms+name) |
| `Fig1_PAGA_branch_hub.py` | PAGA branch-hub / betweenness panel | high (hash-name) |
| `gap_check.py` | continuity check across the ICM/TE-Epi/Hypo axis (no gap) | medium (infer) |

### Figure 2  (`Figure_2/`)

| script | role | confidence |
|---|---|---|
| `marker_dynamics.py` | naive-marker dynamics across stages | high (name) |
| `analysis_progressive_assembly.py` | progressive assembly of the naive programme | medium (name) |
| `figure_LIN28A_MYC_OxPhos.py` | MYC / metabolic panel | low (infer) |

### Figure 3  (`Figure_3/`)

| script | role | confidence |
|---|---|---|
| `analysis_ground_state_module.py` | A,B: ground-state module across the lineage | high (ms+hash-name) |
| `analysis_chen2008_conservation.py` | C,D: Chen 2008 mouse ChIP targets mapped to human orthologues | high (ms+hash-name) |
| `analysis_icm_chromatin_priming.py` | G,H,I: normal-ploidy ICM chromatin accessibility | high (ms+hash-name) |
| `analysis_footprint_network.py` | J,K: TOBIAS footprints (NANOG/GATA6/ESRRB) | high (ms+hash-name) |

### Figure 4  (`Figure_4/`)

| script | role | confidence |
|---|---|---|
| `naive_projection.py` | A,B: cESFW reference-mapping of naive/reset PSC lines | high (ms+name) |
| `Fig_ESRRB_ICMlike_blastoids.py` | C,D,E: ICM% + ESRRB+GATA6+ signature + per-cell ESRRB level | high (name) |
| `Fig4_Dattani_panelFG.py` | F,G: Dattani naive-to-hypoblast timecourse | high (hash-name) |
| `analysis_dattani_icmlike.py` | H: ESRRB-independent ICM Identity Score for the Dattani intermediate | medium (infer) |
| `iyer_blastoid.py` | Iyer blastoid Ctrl vs RapaLink (mTOR inhibition does not induce ESRRB) | high (name) |

### Figure 5  (`Figure_5/`)

| script | role | confidence |
|---|---|---|
| `theory_landscape_test.py` | A-D: empirical quasi-potential, normal-form fit, criticality test | high (ms+hash-name) |
| `theory_bifurcation.py` | bifurcation diagram (fixed points vs pseudotime) | high (ms) |

### Figure 6  (`Figure_6/`)

| script | role | confidence |
|---|---|---|
| `theory_crossspecies_both.py` | A: modality of the ICM order parameter, both order parameters | high (ms+hash-name) |
| `theory_crossspecies.py` | cross-species modality (single order parameter) | high (ms) |
| `mouse_bifurcation.py` | mouse ICM bifurcation | high (ms) |
| `marmoset_bifurcation.py` | marmoset ICM bifurcation | high (ms) |

### Figure 7  (`Figure_7/`)

| script | role | confidence |
|---|---|---|
| `fig7_generalisation.py` | assembles the four-system figure | high (hash-name) |
| `braun_glial_figure.py` | A,B: glial astrocyte/oligodendrocyte decision | high (name) |
| `gut_ISC_transition_test.py` | C,D: intestinal absorptive/secretory decision | high (name) |
| `trophoblast_transition_test.py` | E,F: syncytio/extravillous trophoblast decision | high (hash-name) |
| `lung_injury_test.py` | G,H: AT2->AT1 alveolar conversion (fibrotic lung) | high (name) |

### Figure S1  (`Figure_S1/`)

| script | role | confidence |
|---|---|---|
| `fig_reject_cells.py` | stress/apoptosis/UPR scores + QC metrics by lineage | high (hash-name) |
| `analysis_reject_cells.py` | underlying statistics | medium (name) |

### Figure S2  (`Figure_S2/`)

| script | role | confidence |
|---|---|---|
| `tf_peak_per_stage.py` | TF activity peak per stage | high (name) |

### Figure S3  (`Figure_S3/`)

| script | role | confidence |
|---|---|---|
| `fig_mtor_arms_supp.py` | A,B: ICM->preEPI switch magnitude and mTOR arm decomposition | high (hash-name) |
| `analysis_mtor_arms.py` | underlying arm decomposition | medium (name) |

### Figure S4  (`Figure_S4/`)

| script | role | confidence |
|---|---|---|
| `fig_mouse_E35_E45.py` | mouse E3.5->E4.5 sterol/mTOR delta-scores | high (name) |
| `analysis_mouse_E35_E45_robust.py` | covariate-adjusted robustness | medium (name) |

### Figure S5  (`Figure_S5/`)

| script | role | confidence |
|---|---|---|
| `mouse_sterol_fig.py` | figure build | high (name) |
| `mouse_sterol_fast.py` | paused ESC (GSE81285) + diapause (GSE241462) analysis | high (name) |

### Figure S6  (`Figure_S6/`)

| script | role | confidence |
|---|---|---|
| `figure_nhp_combined.py` | combined NHP cross-species panel | high (name) |
| `figure_nhp_cynomolgus.py` | cynomolgus analysis | high (name) |

### Figure S7  (`Figure_S7/`)

| script | role | confidence |
|---|---|---|
| `robustness_orderparam.py` | three order parameters compared | high (ms+hash-name) |

### Figure S8  (`Figure_S8/`)

| script | role | confidence |
|---|---|---|
| `analysis_model_comparison.py` | three-model BIC comparison | high (ms+hash-name) |

### Figure S9  (`Figure_S9/`)

| script | role | confidence |
|---|---|---|
| `fig_score_distributions.py` | score distributions | high (hash-name) |

### Figure S10  (`Figure_S10/`)

| script | role | confidence |
|---|---|---|
| `fig_landscape_robust.py` | landscape signatures under phi3 vs phi1 | high (ms+hash-name) |
| `theory_landscape_phi3.py` | landscape recomputed under the multigene composite | high (ms) |

### Figure S11  (`Figure_S11/`)

| script | role | confidence |
|---|---|---|
| `fig_dropout_robust.py` | dropout bias test + synthetic benchmark | high (ms+hash-name) |

### Figure S12  (`Figure_S12/`)

| script | role | confidence |
|---|---|---|
| `hema_adult_test.py` | megakaryocyte-erythroid progenitor modality test | high (hash-name) |

### Figure S13  (`Figure_S13/`)

| script | role | confidence |
|---|---|---|
| `pilot_multigene_phi.py` | multigene order parameter per dataset | high (ms) |
| `pilot_variance_decomp.py` | variance decomposition vs technical covariates | high (ms) |
| `analysis_embryo_decomp.py` | embryo-of-origin decomposition | high (ms) |

### Figure S14  (`Figure_S14/`)

| script | role | confidence |
|---|---|---|
| `figure_negative_controls.py` | synthetic resolved-mixture positive control | medium (infer) |

### Figure S15  (`Figure_S15/`)

| script | role | confidence |
|---|---|---|
| `lung_at2at1_test.py` | homeostatic lung (Sikkema 2023) | high (name) |
| `lung_injury_test.py` | fibrotic lung (Habermann 2020) | high (name) |

### Figure S16  (`Figure_S16/`)

| script | role | confidence |
|---|---|---|
| `analysis_footprint_network.py` | TOBIAS footprint network | high (ms+hash-name) |

### Figure S17  (`Figure_S17/`)

| script | role | confidence |
|---|---|---|
| `analysis_icm_vs_naive_atac.py` | ICM vs naive ATAC comparison | high (ms+hash-name) |
| `tobias_fork.py` | TOBIAS macOS fork-mode wrapper | high (ms) |

### Figure S18  (`Figure_S18/`)

| script | role | confidence |
|---|---|---|
| `analysis_ignition_circuit.py` | A-D: ignition criteria, co-expression, wiring, footprints | high (ms+hash-name) |

### Figure S19  (`Figure_S19/`)

| script | role | confidence |
|---|---|---|
| `transition_replicate.py` | per-study replication of programme shapes | medium (name) |
| `tf_replicate.py` | per-study TF-activity replication | medium (name) |
| `epi_replicate.py` | per-study stage-sequence replication | medium (name) |

### Figure S20  (`Figure_S20/`)

| script | role | confidence |
|---|---|---|
| `make_figS20.py` | A-E: per-stage specificity, sens/spec/precision, AUC ranking, panel heatmap | high (authored-here) |
| `icm_signature_specificity.py` | underlying specificity statistics | high (authored-here) |
| `icm_signature_gating.py` | candidate gate search | high (authored-here) |
| `icm_staining_markers.py` | unbiased AUC marker ranking | high (authored-here) |
| `Fig_TableS2_Stirparo_pooling.py` | the OTHER Figure S20 (Stirparo pooling artifact) | high (hash-name) |

### Figure S21  (`Figure_S21/`)

| script | role | confidence |
|---|---|---|
| `make_figS21.py` | A-D: crossover, balance, ignition/transience scatter, per-cell modality | high (authored-here) |
| `ignition_search.py` | genome-wide ignition/transience scan | high (authored-here) |

---

## Known gaps

See `UNRESOLVED.md`. In summary: a small number of *composite* panel images
(figures assembled by pasting several PNGs together) were built by ad-hoc code that
could not be located, and one supplementary figure number is duplicated in the
manuscript. Both are listed explicitly rather than papered over.
