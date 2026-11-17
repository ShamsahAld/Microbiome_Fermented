"""Keep ASVs present in enough cohort samples and write them out for taxonomy.

Writes:
  data/processed/asv_table.biom   cohort samples x prevalent ASVs
  data/processed/asv_seqs.fasta   one sequence per prevalent ASV
"""
from pathlib import Path

import biom
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())

groups = pd.read_csv("data/processed/agp_groups.tsv", sep="\t")
table = biom.load_table("data/raw/agp/agp_16s.biom")
table = table.filter(set(groups["sample_id"]), axis="sample", inplace=False)
print(f"samples: {table.shape[1]}, ASVs: {table.shape[0]}")

min_samples = cfg["asv"]["min_prevalence"] * table.shape[1]
table = table.filter(
    lambda values, id_, md: (values > 0).sum() >= min_samples,
    axis="observation", inplace=False)
print(f"ASVs in >= {cfg['asv']['min_prevalence']:.0%} of samples: {table.shape[0]}")

ids = table.ids(axis="observation")
assert all(set(i.upper()) <= set("ACGTN") for i in ids[:100]), "ASV IDs are not sequences"

Path("data/processed").mkdir(parents=True, exist_ok=True)
with biom.util.biom_open("data/processed/asv_table.biom", "w") as f:
    table.to_hdf5(f, "export_asvs.py")
with open("data/processed/asv_seqs.fasta", "w") as f:
    for i, seq in enumerate(ids):
        f.write(f">{seq}\n{seq}\n")
print("wrote data/processed/asv_table.biom and asv_seqs.fasta")
