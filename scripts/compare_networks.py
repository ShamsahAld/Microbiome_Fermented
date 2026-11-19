"""Compare the consumer and non-consumer networks directly.

Writes results/comparison/summary.tsv and
results/comparison/differential_edges.tsv
"""
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.special import comb

out = Path("results/comparison")
out.mkdir(parents=True, exist_ok=True)
groups = ["consumer", "nonconsumer"]

graphs = {g: nx.read_graphml(f"results/networks/{g}.graphml") for g in groups}
corr = {g: pd.read_csv(f"data/processed/fastspar/{g}/correlation.tsv", sep="\t", index_col=0)
        for g in groups}
nodes = {g: pd.read_csv(f"results/topology/{g}_nodes.tsv", sep="\t").set_index("genus")
         for g in groups}

def edge_set(g):
    return {tuple(sorted(e)) for e in g.edges()}

edges = {g: edge_set(graphs[g]) for g in groups}
shared = edges["consumer"] & edges["nonconsumer"]
union = edges["consumer"] | edges["nonconsumer"]

# Same sign in both networks?
same_sign = sum(
    np.sign(graphs["consumer"].edges[e]["weight"]) == np.sign(graphs["nonconsumer"].edges[e]["weight"])
    for e in shared)

# Hub similarity: degree of every genus in both networks (0 if unconnected)
all_genera = sorted(graphs["consumer"].nodes)
deg_c = [graphs["consumer"].degree(n) for n in all_genera]
deg_n = [graphs["nonconsumer"].degree(n) for n in all_genera]
rho, rho_p = spearmanr(deg_c, deg_n)

# Module similarity: adjusted Rand index on genera connected in both
common = nodes["consumer"].index.intersection(nodes["nonconsumer"].index)
table = pd.crosstab(nodes["consumer"].loc[common, "module"],
                    nodes["nonconsumer"].loc[common, "module"])
sum_cells = comb(table.values, 2).sum()
sum_rows = comb(table.sum(axis=1).values, 2).sum()
sum_cols = comb(table.sum(axis=0).values, 2).sum()
expected = sum_rows * sum_cols / comb(len(common), 2)
ari = (sum_cells - expected) / (0.5 * (sum_rows + sum_cols) - expected)

summary = pd.DataFrame([
    ("edges_consumer", len(edges["consumer"])),
    ("edges_nonconsumer", len(edges["nonconsumer"])),
    ("edges_shared", len(shared)),
    ("edge_jaccard", round(len(shared) / len(union), 3)),
    ("shared_edges_same_sign", f"{same_sign}/{len(shared)}"),
    ("degree_spearman_rho", round(rho, 3)),
    ("degree_spearman_p", f"{rho_p:.1e}"),
    ("genera_connected_in_both", len(common)),
    ("module_adjusted_rand", round(ari, 3)),
], columns=["metric", "value"])
summary.to_csv(out / "summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))

# Differential edges: present in one network only, with r in both groups
rows = []
for g, other in [("consumer", "nonconsumer"), ("nonconsumer", "consumer")]:
    for a, b in sorted(edges[g] - edges[other]):
        rows.append({"only_in": g, "genus_a": a, "genus_b": b,
                     "r_consumer": corr["consumer"].loc[a, b],
                     "r_nonconsumer": corr["nonconsumer"].loc[a, b]})
diff = pd.DataFrame(rows)
diff["r_difference"] = diff["r_consumer"] - diff["r_nonconsumer"]
diff = diff.reindex(diff["r_difference"].abs().sort_values(ascending=False).index)
diff.to_csv(out / "differential_edges.tsv", sep="\t", index=False)

print(f"\nEdges only in consumers: {(diff.only_in == 'consumer').sum()}, "
      f"only in non-consumers: {(diff.only_in == 'nonconsumer').sum()}")
print("\nLargest differences:")
print(diff.head(15).round(3).to_string(index=False))
