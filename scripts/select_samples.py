"""Select AGP stool samples in the consumer / non-consumer groups.

Reads the raw AGP table and metadata, applies the filters in
config/config.yaml, and writes:
  metadata/agp_cohort_ids.txt   one sample ID per line
  data/processed/agp_groups.tsv sample_id <tab> group
"""
from pathlib import Path

import biom
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())["agp"]

table = biom.load_table("data/raw/agp/agp_16s.biom")
depth = pd.Series(table.sum(axis="sample"), index=table.ids(axis="sample"))

meta = pd.read_csv("data/raw/agp/agp_metadata.tsv", sep="\t",
                   index_col=0, dtype=str, low_memory=False)

# Ambiguous samples carry an extra ".suffix" in the metadata;
# map each metadata ID to the ID used in the table.
def to_table_id(sid):
    if sid in depth.index:
        return sid
    base = sid.rsplit(".", 1)[0]
    return base if base in depth.index else None

meta["table_id"] = [to_table_id(s) for s in meta.index]
print(f"metadata rows: {len(meta)}, matched to table: {meta['table_id'].notna().sum()}")
meta = meta.dropna(subset=["table_id"]).drop_duplicates("table_id")

levels = {lvl: "consumer" for lvl in cfg["consumer_levels"]}
levels.update({lvl: "nonconsumer" for lvl in cfg["nonconsumer_levels"]})

keep = meta[meta["sample_type"] == cfg["sample_type"]].copy()
print(f"stool: {len(keep)}")
keep["group"] = keep["fermented_plant_frequency"].map(levels)
keep = keep.dropna(subset=["group"])
print(f"in a group: {len(keep)}")
keep["reads"] = depth.reindex(keep["table_id"]).to_numpy()
keep = keep[keep["reads"] >= cfg["min_reads"]]
print(f">= {cfg['min_reads']} reads: {len(keep)}")
print(keep["group"].value_counts().to_string())

# One sample per person: keep each participant's sample with the most reads
n_missing = keep["host_subject_id"].isna().sum()
if n_missing:
    print(f"dropping {n_missing} samples with no host_subject_id")
keep = keep.dropna(subset=["host_subject_id"])
keep = keep.sort_values("reads", ascending=False).drop_duplicates("host_subject_id")
print(f"one per person: {len(keep)}")

Path("data/processed").mkdir(parents=True, exist_ok=True)
keep[["table_id", "group"]].rename(columns={"table_id": "sample_id"}) \
    .to_csv("data/processed/agp_groups.tsv", sep="\t", index=False)
Path("metadata/agp_cohort_ids.txt").write_text(
    "\n".join(sorted(keep["table_id"])) + "\n")
