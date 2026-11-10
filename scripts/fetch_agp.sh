#!/usr/bin/env bash
# Fetch American Gut Project 16S data (Qiita study 10317) with redbiom.
# Run from the repo root:  pixi run fetch-agp
set -euo pipefail

CTX="Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b"   # pinned redbiom context
IDS="metadata/agp_sample_ids.txt"                    # committed sample list
OUT="data/raw/agp"
CHECKSUMS="metadata/agp_checksums.md5"               # committed md5s

mkdir -p "$OUT"

# Workaround: nltk refuses hard-linked files, and Pixi installs with hard links.
find "$CONDA_PREFIX" -path "*redbiom/assets/nltk_data*" -type f -links +1 \
  -exec sh -c 'cp "$1" "$1.tmp" && mv "$1.tmp" "$1"' _ {} \;

# Skip the download if both files already exist
if [[ ! -s "$OUT/agp_16s.biom" || ! -s "$OUT/agp_metadata.tsv" ]]; then
  echo "Fetching AGP data (this can take a few hours)..."
  rm -f "$OUT"/agp_16s.biom* "$OUT"/agp_metadata.tsv*
  redbiom fetch samples --context "$CTX" --from "$IDS" \
    --resolve-ambiguities most-reads --output "$OUT/agp_16s.biom"
  redbiom fetch sample-metadata --context "$CTX" --from "$IDS" \
    --output "$OUT/agp_metadata.tsv"
fi

# Verify the files match the ones the analysis was built on
md5sum -c "$CHECKSUMS"

# Make raw data read-only
chmod a-w "$OUT"/*
