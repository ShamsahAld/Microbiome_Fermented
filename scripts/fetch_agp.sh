#!/usr/bin/env bash
# Fetch the American Gut Project 16S cohort (Qiita study 10317) with redbiom.
# Run from the repo root:  pixi run fetch-agp
set -euo pipefail

CTX="Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b"   # pinned redbiom context
IDS="metadata/agp_cohort_ids.txt"                    # cohort sample list
OUT="data/raw/agp"
CHECKSUMS="metadata/agp_checksums.md5"               # md5 of metadata + ambiguity files
BIOM_MD5="metadata/agp_biom_data.md5"                # md5 of the table's data only

mkdir -p "$OUT"

# Workaround: nltk refuses hard-linked files, and Pixi installs with hard links.
find "$CONDA_PREFIX" -path "*redbiom/assets/nltk_data*" -type f -links +1 \
  -exec sh -c 'cp "$1" "$1.tmp" && mv "$1.tmp" "$1"' _ {} \;

if [[ ! -s "$OUT/agp_16s.biom" || ! -s "$OUT/agp_metadata.tsv" ]]; then
  echo "Fetching AGP cohort (can take a while)..."
  rm -f "$OUT"/agp_16s.biom* "$OUT"/agp_metadata.tsv*
  redbiom fetch samples --context "$CTX" --from "$IDS" \
    --resolve-ambiguities most-reads --output "$OUT/agp_16s.biom"
  redbiom fetch sample-metadata --context "$CTX" --from "$IDS" \
    --output "$OUT/agp_metadata.tsv"
  # Fixed column/row order so the metadata md5 is the same on every fetch
  python scripts/normalize_metadata.py "$OUT/agp_metadata.tsv"
fi

md5sum -c "$CHECKSUMS"
python scripts/check_biom.py "$OUT/agp_16s.biom" "$BIOM_MD5"

chmod a-w "$OUT"/*
