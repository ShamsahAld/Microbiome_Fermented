"""Turn FastSpar correlations + p-values into one network per group.

An edge joins two genera if |r| >= min_abs_corr and p <= max_pvalue
(both from config/config.yaml). Writes, per group:
  results/networks/{group}.graphml   for analysis
  results/networks/{group}_edges.tsv for reading
and a summary table results/networks/summary.tsv
"""
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())["sparcc"]
out = Path("results/networks")
out.mkdir(parents=True, exist_ok=True)

rows = []
for group in ["consumer", "nonconsumer"]:
    base = Path(f"data/processed/fastspar/{group}")
    r = pd.read_csv(base / "correlation.tsv", sep="\t", index_col=0)
    p = pd.read_csv(base / "pvalues.tsv", sep="\t", index_col=0).loc[r.index, r.columns]

    upper = np.triu(np.ones(r.shape, dtype=bool), k=1)  # each pair once
    keep = upper & (r.abs().values >= cfg["min_abs_corr"]) & (p.values <= cfg["max_pvalue"])
    i, j = np.where(keep)
    edges = pd.DataFrame({
        "source": r.index[i], "target": r.columns[j],
        "r": r.values[i, j], "p": p.values[i, j],
    })
    edges["sign"] = np.where(edges["r"] > 0, "positive", "negative")
    edges.to_csv(out / f"{group}_edges.tsv", sep="\t", index=False)

    g = nx.Graph()
    g.add_nodes_from(r.index)  # keep all genera, even unconnected ones
    for e in edges.itertuples():
        g.add_edge(e.source, e.target, weight=float(e.r), sign=e.sign)
    nx.write_graphml(g, out / f"{group}.graphml")

    connected = [n for n in g if g.degree(n) > 0]
    rows.append({
        "group": group,
        "genera": g.number_of_nodes(),
        "connected_genera": len(connected),
        "edges": g.number_of_edges(),
        "positive": int((edges["sign"] == "positive").sum()),
        "negative": int((edges["sign"] == "negative").sum()),
        "density": round(nx.density(g), 4),
        "mean_degree": round(2 * g.number_of_edges() / g.number_of_nodes(), 2),
        "components": nx.number_connected_components(g.subgraph(connected)),
    })

summary = pd.DataFrame(rows)
summary.to_csv(out / "summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
