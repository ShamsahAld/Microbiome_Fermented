#!/usr/bin/env bash
# Run SparCC (FastSpar) with bootstrap p-values for one group.
# Usage (from the repo root):  bash scripts/run_fastspar.sh consumer
set -euo pipefail

GROUP="$1"
IN="data/processed/fastspar/${GROUP}_counts.tsv"
OUT="data/processed/fastspar/${GROUP}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Settings from config/config.yaml
read -r ITER BOOT SEED < <(python -c "
import yaml
c = yaml.safe_load(open('config/config.yaml'))
print(c['sparcc']['iterations'], c['sparcc']['bootstraps'], c['seed'])")

mkdir -p "$OUT/boot_counts" "$OUT/boot_corr"

echo "[$GROUP] correlations"
fastspar --otu_table "$IN" --correlation "$OUT/correlation.tsv" \
  --covariance "$OUT/covariance.tsv" --iterations "$ITER" \
  --threads "$THREADS" --seed "$SEED" > /dev/null

echo "[$GROUP] $BOOT bootstrap tables"
fastspar_bootstrap --otu_table "$IN" --number "$BOOT" \
  --prefix "$OUT/boot_counts/boot" --seed "$SEED"

echo "[$GROUP] bootstrap correlations"
for f in "$OUT"/boot_counts/boot_*.tsv; do
  name=$(basename "$f")
  fastspar --otu_table "$f" --correlation "$OUT/boot_corr/cor_$name" \
    --covariance "$OUT/boot_corr/cov_$name" --iterations 5 \
    --threads "$THREADS" --seed "$SEED" > /dev/null
done

echo "[$GROUP] p-values"
fastspar_pvalues --otu_table "$IN" --correlation "$OUT/correlation.tsv" \
  --prefix "$OUT/boot_corr/cor_boot_" --permutations "$BOOT" \
  --outfile "$OUT/pvalues.tsv" --threads "$THREADS"

echo "[$GROUP] done: $OUT/correlation.tsv and $OUT/pvalues.tsv"
