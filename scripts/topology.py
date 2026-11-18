"""Hubs, modules and global topology for each group's network.

Modules are found with Louvain on positive edges only (genera that
co-occur), using the seed from config/config.yaml.
Writes, per group, results/topology/{group}_nodes.tsv, plus
results/topology/global.tsv and results/topology/hubs.tsv
"""
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
out = Path("results/topology")
out.mkdir(parents=True, exist_ok=True)

global_rows, hub_rows = [], []
for group in ["consumer", "nonconsumer"]:
    g = nx.read_graphml(f"results/networks/{group}.graphml")
    g = g.subgraph([n for n in g if g.degree(n) > 0]).copy()  # connected genera only

    pos = nx.Graph()
    pos.add_nodes_from(g)
    pos.add_edges_from((u, v, d) for u, v, d in g.edges(data=True) if d["weight"] > 0)
    modules = nx.community.louvain_communities(pos, weight="weight", seed=cfg["seed"])
    module_of = {n: i for i, m in enumerate(sorted(modules, key=len, reverse=True)) for n in m}

    nodes = pd.DataFrame({
        "genus": list(g),
        "degree": [g.degree(n) for n in g],
        "betweenness": pd.Series(nx.betweenness_centrality(g)).reindex(list(g)).values,
        "module": [module_of[n] for n in g],
    }).sort_values(["degree", "betweenness"], ascending=False)
    nodes.to_csv(out / f"{group}_nodes.tsv", sep="\t", index=False)

    top = nodes.head(cfg["topology"]["top_n_hubs"]).assign(group=group)
    hub_rows.append(top)

    largest = max(nx.connected_components(g), key=len)
    global_rows.append({
        "group": group,
        "connected_genera": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "largest_cluster": len(largest),
        "modules": len(modules),
        "modules_3plus": sum(len(m) >= 3 for m in modules),
        "modularity": round(nx.community.modularity(pos, modules, weight="weight"), 3),
        "clustering": round(nx.average_clustering(g), 3),
    })

glob = pd.DataFrame(global_rows)
glob.to_csv(out / "global.tsv", sep="\t", index=False)
hubs = pd.concat(hub_rows)
hubs.to_csv(out / "hubs.tsv", sep="\t", index=False)

print(glob.to_string(index=False))
for group in ["consumer", "nonconsumer"]:
    print(f"\nTop hubs ({group}):")
    print(hubs[hubs.group == group][["genus", "degree", "betweenness", "module"]]
          .round(3).to_string(index=False))
