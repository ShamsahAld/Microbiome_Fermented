"""Split the genus table into one FastSpar input table per group.

Non-consumers are randomly downsampled to the consumer group size (fixed
seed), because network properties depend on sample size.

Writes data/processed/fastspar/{group}_counts.tsv
(FastSpar format: rows = genera, columns = samples, integer counts).
"""
from pathlib import Path

import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())

# Sample IDs look numeric (10317.000001234), so always read them as text
counts = pd.read_csv("data/processed/genus_table.tsv", sep="\t",
                     dtype={"sample_id": str}).set_index("sample_id")
groups = pd.read_csv("data/processed/agp_groups.tsv", sep="\t",
                     dtype={"sample_id": str}).set_index("sample_id")["group"]

missing = set(counts.index) - set(groups.index)
assert not missing, f"{len(missing)} samples in the genus table have no group"
groups = groups.reindex(counts.index)

consumers = counts[groups == "consumer"]
nonconsumers = counts[groups == "nonconsumer"]
if cfg["balance"]["subsample_nonconsumers"]:
    nonconsumers = nonconsumers.sample(n=len(consumers), random_state=cfg["seed"])

out = Path("data/processed/fastspar")
out.mkdir(parents=True, exist_ok=True)
for name, df in [("consumer", consumers), ("nonconsumer", nonconsumers)]:
    table = df.T.astype(int)
    table.index.name = "#OTU ID"
    table.to_csv(out / f"{name}_counts.tsv", sep="\t")
    print(f"{name}: {table.shape[1]} samples x {table.shape[0]} genera")
