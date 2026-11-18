"""Robustness of the consumer vs non-consumer network comparison.

1. Subsampling: rebuild the non-consumer network from several random
   subsets (same size as the consumer group) to see how much it varies.
2. Permutation test: pool consumers with one non-consumer subset, shuffle
   the group labels, rebuild both networks, and record the differences.
   p = share of shuffles with a difference at least as large as the real one.

Writes results/robustness/{subsample,permutation,permutation_summary}.tsv
"""
import os
import subprocess
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
sp, rb, seed = cfg["sparcc"], cfg["robustness"], cfg["seed"]
threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")
out = Path("results/robustness")
out.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv("data/processed/genus_table.tsv", sep="\t",
                     dtype={"sample_id": str}).set_index("sample_id")
groups = pd.read_csv("data/processed/agp_groups.tsv", sep="\t",
                     dtype={"sample_id": str}).set_index("sample_id")["group"]
groups = groups.reindex(counts.index)
consumers = counts[groups == "consumer"]
nonconsumers = counts[groups == "nonconsumer"]
n = len(consumers)


def run(cmd):
    subprocess.run([str(c) for c in cmd], check=True, stdout=subprocess.DEVNULL)


def sparcc_network(df):
    """SparCC with bootstrap p-values; returns the thresholded network."""
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        tmp = Path(tmp)
        table = df.T.astype(int)
        table.index.name = "#OTU ID"
        inp = tmp / "counts.tsv"
        table.to_csv(inp, sep="\t")
        (tmp / "bc").mkdir()
        (tmp / "br").mkdir()
        run(["fastspar", "--otu_table", inp, "--correlation", tmp / "cor.tsv",
             "--covariance", tmp / "cov.tsv", "--iterations", sp["iterations"],
             "--threads", threads, "--seed", seed])
        run(["fastspar_bootstrap", "--otu_table", inp, "--number", sp["bootstraps"],
             "--prefix", tmp / "bc" / "boot", "--seed", seed])
        for f in sorted((tmp / "bc").glob("boot_*.tsv")):
            run(["fastspar", "--otu_table", f, "--correlation", tmp / "br" / f"cor_{f.name}",
                 "--covariance", tmp / "br" / f"cov_{f.name}", "--iterations", 5,
                 "--threads", threads, "--seed", seed])
        run(["fastspar_pvalues", "--otu_table", inp, "--correlation", tmp / "cor.tsv",
             "--prefix", tmp / "br" / "cor_boot_", "--permutations", sp["bootstraps"],
             "--outfile", tmp / "p.tsv", "--threads", threads])
        r = pd.read_csv(tmp / "cor.tsv", sep="\t", index_col=0)
        p = pd.read_csv(tmp / "p.tsv", sep="\t", index_col=0).loc[r.index, r.columns]

    upper = np.triu(np.ones(r.shape, dtype=bool), k=1)
    keep = upper & (r.abs().values >= sp["min_abs_corr"]) & (p.values <= sp["max_pvalue"])
    g = nx.Graph()
    for a, b in zip(*np.where(keep)):
        g.add_edge(r.index[a], r.columns[b], weight=float(r.values[a, b]))
    return g


def metrics(g):
    pos = nx.Graph()
    pos.add_nodes_from(g)
    pos.add_edges_from((u, v, d) for u, v, d in g.edges(data=True) if d["weight"] > 0)
    modules = nx.community.louvain_communities(pos, weight="weight", seed=seed)
    return {
        "connected_genera": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "negative": sum(d["weight"] < 0 for _, _, d in g.edges(data=True)),
        "largest_cluster": len(max(nx.connected_components(g), key=len)),
        "modularity": nx.community.modularity(pos, modules, weight="weight"),
        "clustering": nx.average_clustering(g),
    }


# 1. Observed networks (same subset as the main pipeline: seed + 0)
obs_c = metrics(sparcc_network(consumers))
nonsub = nonconsumers.sample(n=n, random_state=seed)
obs_n = metrics(sparcc_network(nonsub))
obs_diff = {k: obs_c[k] - obs_n[k] for k in obs_c}
print("observed consumer - nonconsumer:", {k: round(v, 3) for k, v in obs_diff.items()}, flush=True)

# 2. Subsampling: how much does the non-consumer network vary?
rows = [dict(repeat=0, **obs_n)]
for k in range(1, rb["subsample_repeats"]):
    sub = nonconsumers.sample(n=n, random_state=seed + k)
    rows.append(dict(repeat=k, **metrics(sparcc_network(sub))))
    print(f"subsample {k + 1}/{rb['subsample_repeats']} done", flush=True)
sub_df = pd.DataFrame(rows)
sub_df.to_csv(out / "subsample.tsv", sep="\t", index=False)

# 3. Permutation test
pool = pd.concat([consumers, nonsub])
rng = np.random.default_rng(seed)
perm_rows = []
for k in range(rb["permutations"]):
    idx = rng.permutation(len(pool))
    ma = metrics(sparcc_network(pool.iloc[idx[:n]]))
    mb = metrics(sparcc_network(pool.iloc[idx[n:]]))
    perm_rows.append({m: ma[m] - mb[m] for m in ma})
    print(f"permutation {k + 1}/{rb['permutations']} done", flush=True)
perm = pd.DataFrame(perm_rows)
perm.to_csv(out / "permutation.tsv", sep="\t", index=False)

summary = pd.DataFrame({
    "metric": list(obs_diff),
    "consumer": [obs_c[m] for m in obs_diff],
    "nonconsumer": [obs_n[m] for m in obs_diff],
    "observed_diff": [obs_diff[m] for m in obs_diff],
    "p_value": [(1 + (perm[m].abs() >= abs(obs_diff[m])).sum()) / (1 + len(perm))
                for m in obs_diff],
})
summary.to_csv(out / "permutation_summary.tsv", sep="\t", index=False)
print("\nSubsampled non-consumer networks:")
print(sub_df.round(3).to_string(index=False))
print("\nPermutation test:")
print(summary.round(3).to_string(index=False))
