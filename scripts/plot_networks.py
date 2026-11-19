"""Draw the consumer and non-consumer networks side by side.

Layout: computed once from positive edges of both networks combined,
so each genus sits in the same place in both panels. Negative edges are
drawn on top (dashed red) but do not affect the layout.
Colour = bacterial family (same in both panels); size = degree.
Writes results/figures/networks.png and .pdf
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd
import yaml

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
out = Path("results/figures")
out.mkdir(parents=True, exist_ok=True)
groups = {"consumer": "Fermented-plant consumers", "nonconsumer": "Non-consumers"}

# Family of each genus label (same labelling rule as collapse_genus.py)
tax = pd.read_csv("data/processed/taxonomy/taxonomy.tsv", sep="\t")
def ranks(t):
    return {r.strip()[:3]: r.strip()[3:] for r in t.split(";") if len(r.strip()) > 3}
def genus_label(n):
    if n.get("g__"):
        return n["g__"]
    for p in ("f__", "o__", "c__", "p__", "d__"):
        if n.get(p):
            return f"{n[p]}_unclassified"
    return "Unassigned"
parsed = tax["Taxon"].map(ranks)
fam = pd.DataFrame({"genus": parsed.map(genus_label),
                    "family": parsed.map(lambda n: n.get("f__") or "Unknown")})
family_of = fam.groupby("genus")["family"].agg(lambda s: s.mode()[0])

graphs = {}
for g in groups:
    full = nx.read_graphml(f"results/networks/{g}.graphml")
    graphs[g] = full.subgraph([n for n in full if full.degree(n) > 0]).copy()

# Shared layout from positive edges only
layout = nx.Graph()
for G in graphs.values():
    layout.add_nodes_from(G)
    layout.add_edges_from(G.edges())
core = max(nx.connected_components(layout), key=len)
layout = layout.subgraph(core).copy()
graphs = {g: G.subgraph([n for n in G if n in core]).copy() for g, G in graphs.items()}
pos = nx.spring_layout(layout, seed=cfg["seed"], weight=None, iterations=500)

# Colour the 8 most common families; group the rest as "Other"
plotted = sorted(layout.nodes)
top = pd.Series([family_of.get(n, "Unknown") for n in plotted]).value_counts().index[:8]
palette = plt.get_cmap("tab10")
fam_colour = {f: palette(i) for i, f in enumerate(top)}
def colour(n):
    return fam_colour.get(family_of.get(n, "Unknown"), "#c7c7c7")

fig, axes = plt.subplots(1, 2, figsize=(18, 9))
for ax, (g, title) in zip(axes, groups.items()):
    G = graphs[g]
    pos_e = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] > 0]
    neg_e = [(u, v) for u, v, d in G.edges(data=True) if d["weight"] < 0]
    nx.draw_networkx_edges(G, pos, edgelist=pos_e, ax=ax,
                           edge_color="#8c8c8c", width=0.7, alpha=0.5)
    nx.draw_networkx_edges(G, pos, edgelist=neg_e, ax=ax, style="dashed",
                           edge_color="#d62728", width=0.8, alpha=0.6)
    deg = dict(G.degree())
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[20 + 22 * deg[n] for n in G],
                           node_color=[colour(n) for n in G],
                           edgecolors="white", linewidths=0.6)
    for n in sorted(deg, key=deg.get, reverse=True)[:6]:
        x, y = pos[n]
        ax.text(x, y + 0.035, n.split("_")[0], fontsize=8, ha="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
    n_neg = len(neg_e)
    ax.set_title(f"{title}\n{G.number_of_nodes()} genera, {G.number_of_edges()} edges "
                 f"({n_neg} negative)", fontsize=12)
    ax.axis("off")

handles = [Line2D([], [], marker="o", ls="", color=c, markersize=8, label=f)
           for f, c in fam_colour.items()]
handles += [Line2D([], [], marker="o", ls="", color="#c7c7c7", markersize=8, label="Other"),
            Line2D([], [], color="#8c8c8c", label="positive"),
            Line2D([], [], color="#d62728", ls="--", label="negative")]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=9)
fig.suptitle("SparCC co-occurrence networks (|r| ≥ 0.3, p ≤ 0.05). "
             "Main connected network; same layout in both panels; colour = family, size = degree", fontsize=11)
fig.tight_layout(rect=(0, 0.06, 1, 0.97))
for ext in ("png", "pdf"):
    fig.savefig(out / f"networks.{ext}", dpi=200, bbox_inches="tight")
print("wrote results/figures/networks.png and .pdf")
