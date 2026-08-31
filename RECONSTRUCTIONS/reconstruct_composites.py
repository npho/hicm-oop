#!/usr/bin/env python
# =============================================================================
#  *** RECONSTRUCTION — THIS IS NOT THE ORIGINAL CODE ***
# =============================================================================
#  Six figures in the manuscript are COMPOSITES: several separately-generated
#  panel images pasted into one. The scripts that produced the individual panels
#  are all present in this package and attributed. The scripts that performed the
#  *compositing* could not be located on disk.
#
#  This file reconstructs that compositing step from the component images. It is
#  written from the evidence (component image names, sizes, aspect ratios and the
#  published captions), NOT recovered from the original source. Treat every
#  output here as "how the panel was most likely assembled", not as provenance.
#
#  It is verifiable: each reconstruction is compared against the published
#  composite still on disk, and the size agreement is reported. A close match
#  supports the reconstruction; a poor match means the real assembly differed
#  (e.g. cropping or manual alignment) and the auditor should rely on the
#  component scripts instead.
#
#  Usage:  python reconstruct_composites.py [--write]
#          (without --write it only reports; with --write it saves
#           *_RECONSTRUCTED.png next to this file)
# =============================================================================
import os
import sys

try:
    from PIL import Image
except ImportError:
    sys.exit('Pillow required:  pip install pillow')

P = '/Users/adelosan/cESFW_analysis/plots'
RR = os.path.join(P, 'revision_results')
OUT = os.path.dirname(os.path.abspath(__file__))
WRITE = '--write' in sys.argv

# figure -> (published composite, [components, top to bottom], note)
JOBS = {
    'Figure 6  Fig6_crossspecies_plus_lineage': (
        None,  # published composite exists only inside the .docx
        [f'{RR}/Fig_crossspecies_both.png',
         f'{RR}/Fig_downstream_resolution.png'],
        'Caption: (A) mouse vs primate ICM modality; (B) resolution along the '
        'human lineage. A = theory_crossspecies_both.py, B = downstream '
        'resolution. Order inferred from the caption.'),
    'Figure S15  Fig_lung_combined': (
        f'{RR}/Fig_lung_combined.png',
        [f'{RR}/Fig_lung_at2at1.png',
         f'{RR}/Fig_lung_injury.png'],
        'Homeostatic (Sikkema 2023) above, fibrotic (Habermann 2020) below. '
        'Component sizes 1950x780 each; published composite 1950x1734, i.e. '
        'two stacked panels plus ~174px of spacing/labels.'),
    'Figure 3  Fig3_merged_withfootprints': (
        f'{P}/Fig3_merged_withfootprints.png',
        [f'{RR}/Fig3_part2_GHIJ.png',
         f'{RR}/Fig_ICM_dualcompetence_footprints.png'],
        'Recorded during manuscript prep: the footprint panel (K) was '
        'composited BELOW the G-J image. Published composite is 2788x1915 vs '
        'G-J alone at 2788x1668, so the footprint strip was scaled/cropped to '
        'about 247px — the reconstruction will therefore NOT match exactly.'),
    'Figure 4  Fig4_merged_ABCDEFG': (
        f'{P}/Fig4_merged_ABCDEFG.png',
        [f'{P}/media_v34/media/image7.png',
         f'{P}/Fig_ESRRB_ICMlike_blastoids.png'],
        'Recorded during manuscript prep: composite of the A/B composition '
        'raster (media_v34 image7.png, with its top ~100px negative header '
        'CROPPED OFF) above the C/D/E figure. The A/B sub-titles are baked into '
        'the raster and its own source script was never located.'),
}

# figures whose "missing" image is a single panel, not a composite
SINGLE = {
    'Figure S13  Fig_cesfw_independent': (
        f'{RR}/Fig_cesfw_independent.png',
        'Not a composite. Produced by the cESFW-independence analysis '
        '(pilot_multigene_phi.py / pilot_variance_decomp.py / '
        'analysis_embryo_decomp.py, all in Figure_S13/). The plotting call was '
        'not located; re-run those analyses to regenerate.'),
    'Figure S14  Fig_negative_control': (
        f'{RR}/Fig_negative_control.png',
        'Not a composite. Attributed to figure_negative_controls.py by '
        'inference only (evidence code: infer). VERIFY before relying on it.'),
}


def stack(paths, pad=40, bg=(255, 255, 255)):
    ims = [Image.open(p).convert('RGB') for p in paths]
    w = max(i.width for i in ims)
    scaled = []
    for i in ims:
        if i.width != w:
            i = i.resize((w, round(i.height * w / i.width)), Image.LANCZOS)
        scaled.append(i)
    h = sum(i.height for i in scaled) + pad * (len(scaled) - 1)
    canvas = Image.new('RGB', (w, h), bg)
    y = 0
    for i in scaled:
        canvas.paste(i, (0, y))
        y += i.height + pad
    return canvas


print('=' * 74)
print('RECONSTRUCTION OF COMPOSITE FIGURES  (not original code)')
print('=' * 74)
for name, (published, comps, note) in JOBS.items():
    print(f'\n--- {name}')
    missing = [c for c in comps if not os.path.exists(c)]
    if missing:
        print('    SKIP, missing components:')
        for m in missing:
            print('      ', m)
        continue
    for c in comps:
        im = Image.open(c)
        print(f'    component  {im.size}  {os.path.relpath(c, P)}')
    rec = stack(comps)
    print(f'    reconstructed -> {rec.size}')
    if published and os.path.exists(published):
        pub = Image.open(published)
        dw = 100 * abs(rec.width - pub.width) / pub.width
        dh = 100 * abs(rec.height - pub.height) / pub.height
        verdict = ('close' if dw < 3 and dh < 8 else
                   'DIFFERS - real assembly involved cropping/scaling')
        print(f'    published     {pub.size}   width diff {dw:.1f}%, '
              f'height diff {dh:.1f}%  -> {verdict}')
    else:
        print('    published composite not on disk (lives only inside the .docx)')
    print(f'    NOTE: {note}')
    if WRITE:
        dest = os.path.join(OUT, name.split()[-1] + '_RECONSTRUCTED.png')
        rec.save(dest, dpi=(300, 300))
        print(f'    wrote {dest}')

print('\n' + '=' * 74)
print('NOT COMPOSITES - single panels whose plotting call was not located')
print('=' * 74)
for name, (img, note) in SINGLE.items():
    ok = 'on disk' if os.path.exists(img) else 'ABSENT'
    print(f'\n--- {name}  ({ok})')
    print(f'    {note}')

print('\nDone. Nothing here should be cited as the original figure code.')
