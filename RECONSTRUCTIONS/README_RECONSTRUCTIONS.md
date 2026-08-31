# RECONSTRUCTIONS — not original code

**Nothing in this folder is recovered source.** These are reconstructions written from
evidence (component image names, pixel dimensions, aspect ratios, the published captions,
and notes recorded during manuscript preparation). They exist so an auditor can see how a
composite figure was *most likely* assembled and check that assumption, not so anyone can
claim this is what was run.

Every figure's **component panels are attributed and present** in the `Figure_*/` folders.
Only the final pasting step is reconstructed.

Run `python reconstruct_composites.py` to print the comparison; add `--write` to save
`*_RECONSTRUCTED.png` for visual comparison against the published panel.

## Verification results

The script compares each reconstruction against the published composite still on disk:

| figure | width diff | height diff | verdict |
|---|---|---|---|
| **Figure S15** `Fig_lung_combined` | 0.0% | 7.7% | **close** — reconstruction is credible |
| Figure 3 `Fig3_merged_withfootprints` | 0.0% | 37.8% | **differs** — footprint strip was scaled/cropped |
| Figure 4 `Fig4_merged_ABCDEFG` | 0.0% | 24.0% | **differs** — real assembly cropped/scaled |
| Figure 6 `Fig6_crossspecies_plus_lineage` | — | — | published composite exists only inside the .docx |

Width matches exactly (0.0%) in all three testable cases, which is good evidence that the
correct component images have been identified. The height mismatches mean the original
assembly cropped or rescaled panels in ways this simple stack does not reproduce.

**For Figure 4 specifically**, the A/B component (`media_v34/media/image7.png`, 4991px wide)
matches the published composite width exactly (4991px), confirming it is the right source.
Its own generating script was never located, and its panel sub-titles ("not ICM", "do not
populate") are baked into the raster.

## What this means for an audit

- For **numerical claims**, ignore this folder. Every number in the paper comes from the
  attributed analysis scripts in `Figure_*/`, which are original.
- For **figure layout**, treat the composites as unreproduced. The honest statement is that
  the panels are reproducible but the final assembly was manual.
- **Figure S13** (`Fig_cesfw_independent`) and **Figure S14** (`Fig_negative_control`) are
  *not* composites — they are single panels whose plotting call was not located. S14 is
  attributed by inference only and should be verified first.
