# Unresolved items for the auditor

Listed explicitly rather than silently omitted.

## 1. Duplicate supplementary figure number (manuscript defect)

The manuscript contains **two figures numbered "Figure S20"**:

- `Figure S20. The Stirparo bimodal ICM is a pooling artifact.` — appears **out of
  sequence**, between S8 and S9 in document order. Produced by
  `Fig_TableS2_Stirparo_pooling.py`.
- `Figure S20. Specificity of the ICM signature, and why TFCP2L1 cannot localise it.`
  — appears after S19, in sequence.

Both are in `Figure_S20/`. This needs renumbering before submission; it is a
manuscript defect, not an attribution failure.

## 2. Composite panel images with no locatable build script

**Now partly addressed** — see `RECONSTRUCTIONS/`, which reconstructs the compositing from the component panels and reports how well each matches the published figure (Figure S15 close at 7.7% height difference; Figures 3 and 4 differ, meaning the original assembly cropped or rescaled). The reconstructions are explicitly NOT original code.


Several figures are composites assembled from multiple PNGs. The component analyses
are all present, but the *compositing* step for these was not found on disk:

| composite image | status |
|---|---|
| `Fig6_crossspecies_plus_lineage.png` | components attributed (theory_crossspecies_both.py etc.); compositing code not found |
| `Fig_lung_combined.png` | components attributed (lung_at2at1_test.py, lung_injury_test.py); compositing code not found |
| `Fig_cesfw_independent.png` | analyses attributed (pilot_*.py, analysis_embryo_decomp.py); compositing code not found |
| `Fig_negative_control.png` | attributed to figure_negative_controls.py by inference only |
| `Fig3_merged_withfootprints.png` | composited in-place by `build_v37.py` (a document builder, not an analysis script) |
| `Fig4_merged_ABCDEFG.png` | composited in-place by `build_v44.py` (document builder) |

Note that `Fig4_merged_ABCDEFG.png` and the Figure 4 A/B sub-panels are partly baked
into a raster whose original source script has not been located — this was already
known and recorded during manuscript preparation.

## 3. Figure 2 and Figure S13 embedded images had no byte-identical match on disk

The embedded images for Figure 2 (two panels) and Figure S13 did not md5-match any PNG
currently on disk, meaning the source PNG has since been regenerated or overwritten.
Scripts for these are attributed by `name`/`infer` and should be re-run to confirm they
reproduce the published panel.

## 4a. CORRECTED — Figure S14 was MISATTRIBUTED (found and fixed 2026-08-28)

The first build of this package attributed Figure S14 to `figure_negative_controls.py`.
**That was wrong.** That script builds naive-PSC projection controls
(t2iLGo/PXGL/Mazid/Bi/EPSC) and writes `negative_controls_extended.png` — unrelated to
S14's Hartigan dip test on Stirparo-labelled embryo lineages. It has been removed.

The real producer **does not exist on disk**. Exhaustive search: no `.py`, `.ipynb`,
`.sh`, `.txt` or `.md` file anywhere under `cESFW_analysis/` or `cesfw_runs/` references
`Fig_negative_control.png` (written 2026-07-19 09:25), and no Python file was modified
that day except the unrelated generalisation tests. The same is true of
`Fig_cesfw_independent.png` (Figure S13, written 26 minutes earlier) — though S13's
underlying analyses ARE attributed, via the Methods.

In place of a false attribution, `Figure_S14/RECONSTRUCTED_negative_control.py`
re-implements the analysis from the caption and Methods. **It reproduces the published
result**: pooled real epiblast+hypoblast is bimodal (dip p=0.014 vs published 0.02) while
ICM (0.856 vs 0.95), epiblast, hypoblast and trophectoderm are each unimodal, with
trophectoderm n=993 matching exactly.

> **One discrepancy worth noting:** the caption states trophectoderm "ΔBIC>1000"; the
> reconstruction gets **963.7**. Close, and it does not affect the argument (the point is
> that ΔBIC is uninformative at large n while the dip test still reads unimodal), but the
> published figure is marginally overstated unless different GMM settings were used.

## 4b. Attribution marked `infer`

These entries were attributed from caption content plus reading the script, without
hash or Methods confirmation. They are the ones to check first:

- **Figure 1** — `gap_check.py` (continuity check across the ICM/TE-Epi/Hypo axis (no gap))
- **Figure 2** — `figure_LIN28A_MYC_OxPhos.py` (MYC / metabolic panel)
- **Figure 4** — `analysis_dattani_icmlike.py` (H: ESRRB-independent ICM Identity Score for the Dattani intermediate)

## 5. Superseded material deliberately excluded

The earlier `paper_code_audit/08_exploratory_or_superseded/` scripts are **not**
included here, because none of their outputs appear in the current manuscript. This
includes the discarded from-scratch linear-PCA trajectory (which flattened a branching
topology onto one axis and was abandoned), the NHP sterol conservation test (did not
conserve cleanly), and the bulk-ATAC ESRRB motif enrichment (non-significant). They
remain available in that folder for transparency.

## 6. Data provenance (added after the first build)

`DATA_PROVENANCE.md` and `fetch_data.sh` were added because the scripts alone are not
sufficient to obtain the data: only 11 accessions appear anywhere in the code, only two
scripts contain a download step, and all paths are absolute to one workstation. Six
datasets remain **portal-only** and cannot be scripted (cESFW atlas on Mendeley, Burclaff
gut on UCSC Cell Browser, Arutyunyan placenta and Sikkema lung on CELLxGENE, Braun brain,
Zeng marrow), plus the Dattani 2024 processed object which came from the authors.

Note also that several PSC inputs (`Guo2021_sc_Counts.csv`, `Polo2020_Naive_QC_counts.csv`,
`Moya2022_QC_counts.csv`, `Bi2022_5iLAF_nESC_QC_counts.csv`, `Mazid2022_sc_Counts.csv`,
`HEPSC_Counts.csv`) are **QC-filtered derivatives**, and the QC step itself is not in this
package.
