#!/usr/bin/env bash
# Download the primary public datasets for the ICM figure scripts.
#
# NOT fully automatic: several datasets sit behind interactive portals
# (Mendeley Data, CELLxGENE, UCSC Cell Browser) and must be fetched by hand.
# Those are listed at the bottom. Everything from GEO/ArrayExpress is scripted.
#
# Usage:  bash fetch_data.sh /path/to/data
set -euo pipefail
DEST="${1:-./data}"
mkdir -p "$DEST"
cd "$DEST"
echo "Downloading into: $(pwd)"

geo () {  # $1 = accession
  local acc="$1"; local stub="${acc:0:${#acc}-3}nnn"
  echo "--- $acc"
  mkdir -p "$acc" && cd "$acc"
  # -r over the suppl dir; drops the FTP path prefix
  wget -nv -r -nH --cut-dirs=5 --no-parent -A '*.gz,*.txt,*.csv,*.tsv,*.bed,*.narrowPeak' \
      "https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${acc}/suppl/" || \
      echo "    (nothing retrieved for $acc - check the accession)"
  cd ..
}

for acc in GSE171820 GSE123046 GSE241462 GSE81285 GSE267302 GSE101571 \
           GSE11431 GSE74767 GSE106481 GSE135893; do
  geo "$acc"
done

echo "--- E-MTAB-7078 (marmoset, Boroviak 2018)"
wget -nv -r -nH --cut-dirs=6 --no-parent -A '*.txt,*.csv,*.tsv,*.gz' \
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/078/E-MTAB-7078/Files/" || true

echo "--- E-MTAB-6819 (Messmer 2019 t2iLGo)"
wget -nv -r -nH --cut-dirs=6 --no-parent -A '*.txt,*.csv,*.tsv,*.gz' \
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/819/E-MTAB-6819/Files/" || true

cat <<'EOF'

=======================================================================
MANUAL DOWNLOADS STILL REQUIRED
=======================================================================
1. cESFW human embryo atlas  (needed by ~43 of 65 scripts)
   https://data.mendeley.com/datasets/34td4ds2r9
   -> place Human_Embryo_Counts.csv, Human_Sample_Info.csv,
      Saved_cESFW_Genes.npy, UMAP_Human_Embryo_Embedding.npy in $DEST

2. Burclaff 2022 intestine
   https://cells.ucsc.edu/?ds=intestine-epithelium
   -> gut_matrix.mtx.gz, gut_features.tsv.gz, gut_barcodes.tsv.gz, gut_meta.tsv

3. Arutyunyan 2023 placenta  and  Sikkema 2023 lung atlas
   https://cellxgene.cziscience.com/   (search by author/year)

4. Braun 2023 first-trimester brain  -> braun_meta.tsv, braun_markers.tsv
5. Zeng 2025 bone marrow reference   (see manuscript Methods)
6. Dattani 2024 processed object     (from the authors)
7. Gao 2018 ATAC peaks               (icm_2pn_peaks.bed.gz)

Gene sets download themselves on first use:
   MSigDB Hallmark 2020 via gseapy  (cached in ~/.cache/gseapy/)
   CollecTRI via  decoupler.get_collectri(organism="human")
=======================================================================
EOF
