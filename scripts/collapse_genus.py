"""Sum ASV counts into genera and keep genera present in enough samples.

Writes data/processed/genus_table.tsv (samples x genera, raw counts).
ASVs without a genus are labelled by their deepest known rank,
e.g. "Lachnospiraceae_unclassified".
"""
from pathlib import Path

import biom
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())

tax = pd.read_csv("data/processed/taxonomy/taxonomy.tsv", sep="\t", index_col=0)

def genus_label(taxon):
    ranks = [r.strip() for r in taxon.split(";")]
    names = {r[:3]: r[3:] for r in ranks if len(r) > 3}
    if names.get("g__"):
        return names["g__"]
    for prefix in ("f__", "o__", "c__", "p__", "d__"):
        if names.get(prefix):
            return f"{names[prefix]}_unclassified"
    return "Unassigned"

tax["genus"] = tax["Taxon"].map(genus_label)

table = biom.load_table("data/processed/asv_table.biom")
counts = pd.DataFrame(table.matrix_data.toarray().T,
                      index=table.ids("sample"),
                      columns=table.ids("observation"))
genus_counts = counts.T.groupby(tax["genus"].reindex(counts.columns)).sum().T
print(f"genera before filtering: {genus_counts.shape[1]}")

min_samples = cfg["genus"]["min_prevalence"] * len(genus_counts)
keep = (genus_counts > 0).sum() >= min_samples
genus_counts = genus_counts.loc[:, keep]
kept_reads = genus_counts.values.sum() / counts.values.sum()
print(f"genera in >= {cfg['genus']['min_prevalence']:.0%} of samples: {genus_counts.shape[1]}")
print(f"reads kept: {kept_reads:.1%}")

genus_counts.to_csv("data/processed/genus_table.tsv", sep="\t")
print("wrote data/processed/genus_table.tsv")
