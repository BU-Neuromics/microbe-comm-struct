#!/usr/bin/env python3
"""
visualize_results.py
Visualization of NetCoMi microbial co-occurrence network results.
Nodes are TAXA; edges represent significant co-occurrence associations.

Produces a PDF with the following figures:
  1.  Title page
  2.  Network overview: edge counts, density, LCC size per method
  3.  Global network properties (clustering coef, modularity, etc.)
  4.  Network graphs — one per method (spring layout, edge sign colored)
  5.  Connected components: size distribution per method
  6.  Cluster/community composition: taxon counts per cluster per method
  7.  Top hub taxa by degree centrality per method
  8.  Top hub taxa by betweenness centrality per method
  9.  Degree distribution per method
  10. Closeness centrality top taxa per method
  11. Positive vs negative edge breakdown per method
  12. Method agreement: edge Jaccard overlap heatmap
  13. Method agreement: degree correlation heatmap
  14. Shared hub taxa across methods (upset-style bar chart)
  15. Key findings summary
"""

import argparse
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import networkx as nx
from scipy import stats

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Visualize NetCoMi microbial co-occurrence network results.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-r",
    "--results-dir",
    default="results",
    metavar="DIR",
    help="Results directory: either a multi-level root (contains class/, family/, …) "
    "or a single-level dir (contains sparcc/, spring/, …).",
)
parser.add_argument(
    "-o",
    "--output",
    default="network_analysis_visualizations.pdf",
    metavar="FILE",
    help="Output PDF file path.",
)
parser.add_argument(
    "--top-n",
    type=int,
    default=30,
    metavar="N",
    help="Number of top taxa to show in centrality bar charts.",
)
args = parser.parse_args()

RESULTS_DIR = Path(args.results_dir)
OUTPUT_PDF = args.output
TOP_N = args.top_n

ALL_TAX_LEVELS = ["class", "order", "family", "genus"]
ALL_METHODS = [
    "sparcc",
    "spearman",
    "cclasso",
    "spieceasi",
    "spring",
    "lupine",
    "ccrepe",
    "propr",
]

METHOD_COLORS = {
    "sparcc": "#E41A1C",
    "spearman": "#FF7F00",
    "cclasso": "#4DAF4A",
    "spieceasi": "#377EB8",
    "spring": "#984EA3",
}
LEVEL_COLORS = {
    "class": "#66C2A5",
    "order": "#FC8D62",
    "family": "#8DA0CB",
    "genus": "#E78AC3",
}
POS_COLOR = "#1565C0"  # positive association edges
NEG_COLOR = "#C62828"  # negative association edges


# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────
def read_safe(path):
    p = Path(path)
    return pd.read_csv(p) if p.exists() else None


def read_matrix(path):
    p = Path(path)
    if not p.exists():
        return None
    df = pd.read_csv(p, index_col=0)
    df.columns = df.columns.str.strip()
    df.index = df.index.str.strip()
    return df.astype(float)


def save_page(pdf, fig):
    if fig is not None:
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def method_color(mth):
    return METHOD_COLORS.get(mth, "#888888")


# ─────────────────────────────────────────────────────────────
# Detect single-level vs multi-level layout
# ─────────────────────────────────────────────────────────────
print("=" * 55)
print("NetCoMi Network Results Visualization")
print("=" * 55)

_subdirs = {p.name for p in RESULTS_DIR.iterdir() if p.is_dir()}
# A dir is "single-level" if it directly contains method subdirs (not tax-level subdirs).
# We treat any subdir that is not a known tax level and not "comparison" as a potential method.
_non_tax_subdirs = _subdirs - set(ALL_TAX_LEVELS) - {"comparison"}
_is_single = bool(_non_tax_subdirs) and not bool(_subdirs & set(ALL_TAX_LEVELS))

if _is_single:
    _inferred = RESULTS_DIR.name if RESULTS_DIR.name in ALL_TAX_LEVELS else "unknown"
    level_dirs = {_inferred: RESULTS_DIR}
    TAX_LEVELS = [_inferred]
    print(f"Single-level mode: '{_inferred}' in {RESULTS_DIR}")
else:
    level_dirs = {
        lvl: RESULTS_DIR / lvl for lvl in ALL_TAX_LEVELS if (RESULTS_DIR / lvl).is_dir()
    }
    TAX_LEVELS = list(level_dirs.keys())
    print(f"Multi-level mode: {TAX_LEVELS}")

# ─────────────────────────────────────────────────────────────
# Load all data
# ─────────────────────────────────────────────────────────────
print("Loading data ...")

data = {}  # data[(lvl, mth)] = dict of DataFrames


def _synthesize_centralities(node_attrs):
    """Build a minimal centralities DataFrame from node_attributes when the
    standard _centralities.csv is absent (e.g. lupine output).

    node_attrs columns: taxon, degree, [mean_abundance, prevalence, ...]
    Returns a DataFrame with columns: taxon, degree, betweenness, closeness
    where betweenness and closeness are NaN (unavailable).
    """
    if node_attrs is None:
        return None
    df = node_attrs[["taxon", "degree"]].copy()
    df["betweenness"] = float("nan")
    df["closeness"] = float("nan")
    return df


for lvl, lvl_dir in level_dirs.items():
    # Discover all method subdirectories present (not just a fixed allowlist)
    candidate_dirs = [
        p for p in lvl_dir.iterdir() if p.is_dir() and p.name not in {"comparison"}
    ]
    for mth_dir in candidate_dirs:
        mth = mth_dir.name
        prefix = mth_dir / f"result_{mth}"
        centralities = read_safe(f"{prefix}_centralities.csv")
        node_attrs = read_safe(f"{prefix}_node_attributes.csv")
        entry = {
            "centralities": centralities,
            "clusters": read_safe(f"{prefix}_clusters.csv"),
            "global_props": read_safe(f"{prefix}_global_properties.csv"),
            "adjacency_matrix": read_matrix(f"{prefix}_adjacency_matrix.csv"),
            "association_matrix": read_matrix(f"{prefix}_association_matrix.csv"),
            # Extra files present for some methods (e.g. lupine)
            "node_attributes": node_attrs,
            "edge_list": read_safe(f"{prefix}_edge_list.csv"),
            "pvalue_matrix": read_matrix(f"{prefix}_pvalue_matrix.csv"),
        }
        # When centralities are absent, synthesize from node_attributes if available
        if entry["centralities"] is None and node_attrs is not None:
            entry["centralities"] = _synthesize_centralities(node_attrs)
        # Store if we have at least an adjacency or association matrix
        if (
            entry["adjacency_matrix"] is not None
            or entry["association_matrix"] is not None
        ):
            data[(lvl, mth)] = entry

# Comparison files
comp_data = {}
for lvl, lvl_dir in level_dirs.items():
    comp_dir = lvl_dir / "comparison"
    comp_data[lvl] = {
        "basic_stats": read_safe(comp_dir / "comparison_basic_stats.csv"),
        "edge_overlap": read_matrix(comp_dir / "comparison_edge_overlap.csv"),
        "deg_corr": read_matrix(comp_dir / "comparison_degree_correlation.csv"),
        "hub_overlap": read_matrix(comp_dir / "comparison_hub_overlap.csv"),
    }

# Methods actually present (with data)
present_methods = sorted(
    {mth for (_, mth) in data.keys()},
    key=lambda m: ALL_METHODS.index(m) if m in ALL_METHODS else 99,
)

print(f"Methods with data: {present_methods}")
print(f"Data loaded: {len(data)} (level, method) combinations.")

# ─────────────────────────────────────────────────────────────
# Global properties helper: parse wide property table
# ─────────────────────────────────────────────────────────────
PROPS = {
    "lccSize1": "LCC size (taxa)",
    "lccSizeRel1": "LCC size (relative)",
    "clustCoef1": "Clustering coefficient",
    "modularity1": "Modularity",
    "natConnect1": "Natural connectivity",
    "avPath1": "Avg path length",
    "density1": "Density",
    "pep1": "Positive edge %",
}


def get_global_prop(lvl, mth, prop):
    entry = data.get((lvl, mth))
    if entry is None:
        return np.nan
    gp = entry["global_props"]
    if gp is None:
        return np.nan
    row = gp[gp["property"] == prop]
    if row.empty:
        return np.nan
    return pd.to_numeric(row["value"].values[0], errors="coerce")


# ─────────────────────────────────────────────────────────────
# Build global props dataframe
# ─────────────────────────────────────────────────────────────
gp_rows = []
for lvl, mth in data:
    row = {"level": lvl, "method": mth}
    for prop in PROPS:
        row[prop] = get_global_prop(lvl, mth, prop)
    gp_rows.append(row)
gp_df = pd.DataFrame(gp_rows)

# ─────────────────────────────────────────────────────────────
# Build combined centralities dataframe
# ─────────────────────────────────────────────────────────────
cent_frames = []
for (lvl, mth), entry in data.items():
    if entry["centralities"] is None:
        continue
    c = entry["centralities"].copy()
    c["level"] = lvl
    c["method"] = mth
    cent_frames.append(c)
cent_df = pd.concat(cent_frames, ignore_index=True) if cent_frames else pd.DataFrame()

# ─────────────────────────────────────────────────────────────
# Open PDF
# ─────────────────────────────────────────────────────────────
print(f"Writing to {OUTPUT_PDF} ...")

with PdfPages(OUTPUT_PDF) as pdf:
    # ── Fig 1: Title page ────────────────────────────────────────────────
    print("  Fig 1: title page")
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis("off")
    level_str = ", ".join(TAX_LEVELS)
    method_str = ", ".join(present_methods)
    ax.text(
        0.5,
        0.75,
        "NetCoMi Microbial Co-occurrence Networks",
        ha="center",
        va="center",
        fontsize=20,
        fontweight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.65,
        "Results Visualization Report",
        ha="center",
        va="center",
        fontsize=15,
        color="#555555",
        transform=ax.transAxes,
    )
    details = (
        f"Taxonomic level(s): {level_str}\n"
        f"Methods: {method_str}\n\n"
        "Nodes = taxa  |  Edges = significant co-occurrence associations\n"
        "Positive edge = co-occurrence  |  Negative edge = mutual exclusion"
    )
    ax.text(
        0.5,
        0.42,
        details,
        ha="center",
        va="center",
        fontsize=11,
        linespacing=2.0,
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.10,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        ha="center",
        va="center",
        fontsize=9,
        color="grey",
        transform=ax.transAxes,
    )
    save_page(pdf, fig)

    # ── Fig 2: Network overview bar charts ──────────────────────────────
    print("  Fig 2: network overview")

    # Gather basic stats
    bs_frames = [
        d["basic_stats"].assign(level=lvl)
        for lvl, d in comp_data.items()
        if d["basic_stats"] is not None
    ]
    if bs_frames:
        bs_df = pd.concat(bs_frames, ignore_index=True)
        bs_df = bs_df[bs_df["method"].isin(present_methods)]

        n_lvls = len(TAX_LEVELS)
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle("Network Overview", fontsize=15, fontweight="bold")

        metrics = [
            ("n_nodes", "Number of taxa (nodes)", False, None),
            ("n_edges", "Number of edges", False, None),
            ("edge_density", "Edge density", True, 1.0),
            ("lcc_size", "Largest connected component (taxa)", False, None),
        ]
        bar_width = 0.7 / max(len(TAX_LEVELS), 1)

        for ax, (col, title, add_hline, hline_val) in zip(axes.flat, metrics):
            for i, lvl in enumerate(TAX_LEVELS):
                sub = bs_df[bs_df["level"] == lvl].copy()
                if sub.empty:
                    continue
                x = [
                    present_methods.index(m) + i * bar_width
                    for m in sub["method"]
                    if m in present_methods
                ]
                y = sub[col].clip(upper=1.2 if col == "edge_density" else None)
                colors = [
                    method_color(m) for m in sub["method"] if m in present_methods
                ]
                ax.bar(
                    x,
                    y,
                    width=bar_width * 0.9,
                    color=colors
                    if len(TAX_LEVELS) == 1
                    else LEVEL_COLORS.get(lvl, "#888"),
                    label=lvl,
                    alpha=0.88,
                )
            if add_hline and hline_val is not None:
                ax.axhline(hline_val, color="grey", linestyle="--", linewidth=0.9)
            cx = np.arange(len(present_methods)) + bar_width * (len(TAX_LEVELS) - 1) / 2
            ax.set_xticks(cx)
            ax.set_xticklabels(present_methods, rotation=25, ha="right", fontsize=9)
            ax.set_title(title, fontsize=10, fontweight="bold")
            if len(TAX_LEVELS) > 1:
                ax.legend(fontsize=8)
            else:
                # Color x-tick labels by method color when single level
                for tick, mth in zip(ax.get_xticklabels(), present_methods):
                    tick.set_color(method_color(mth))

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 3: Global network properties ────────────────────────────────
    print("  Fig 3: global network properties")
    props_panels = [
        ["lccSize1", "lccSizeRel1", "clustCoef1", "modularity1"],
        ["natConnect1", "avPath1", "density1", "pep1"],
    ]
    for panel_idx, panel_props in enumerate(props_panels):
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle(
            f"Global Network Properties — panel {panel_idx + 1}/2",
            fontsize=14,
            fontweight="bold",
        )
        for ax, prop in zip(axes.flat, panel_props):
            label = PROPS[prop]
            sub = gp_df[["level", "method", prop]].dropna(subset=[prop])
            if sub.empty:
                ax.axis("off")
                ax.set_title(label + "\n(no data)", fontsize=9)
                continue
            if len(TAX_LEVELS) == 1:
                lvl = TAX_LEVELS[0]
                s = sub[sub["level"] == lvl].sort_values("method")
                ax.bar(
                    range(len(s)),
                    s[prop],
                    color=[method_color(m) for m in s["method"]],
                    alpha=0.88,
                )
                ax.set_xticks(range(len(s)))
                ax.set_xticklabels(s["method"], rotation=25, ha="right", fontsize=9)
                for tick, mth in zip(ax.get_xticklabels(), s["method"]):
                    tick.set_color(method_color(mth))
            else:
                for mth in present_methods:
                    m_sub = sub[sub["method"] == mth]
                    xs = [
                        TAX_LEVELS.index(l) for l in m_sub["level"] if l in TAX_LEVELS
                    ]
                    ys = [
                        m_sub.loc[m_sub["level"] == l, prop].values[0]
                        for l in m_sub["level"]
                        if l in TAX_LEVELS
                    ]
                    ax.plot(
                        xs,
                        ys,
                        "o-",
                        color=method_color(mth),
                        label=mth,
                        linewidth=1.5,
                        markersize=7,
                    )
                ax.set_xticks(range(len(TAX_LEVELS)))
                ax.set_xticklabels(TAX_LEVELS, rotation=20, ha="right")
                ax.legend(fontsize=7)
            ax.set_title(label, fontsize=10, fontweight="bold")
            ax.set_ylabel("Value")
        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 4: Network graphs ────────────────────────────────────────────
    print("  Fig 4: network graphs")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        for mth in methods_here:
            entry = data[(lvl, mth)]
            assoc = entry["association_matrix"]
            adj = entry["adjacency_matrix"]
            if assoc is None or adj is None:
                continue

            # Mask association by adjacency (keep only retained edges)
            assoc_np = assoc.values.copy()
            adj_np = adj.values.copy()
            np.fill_diagonal(assoc_np, 0)
            np.fill_diagonal(adj_np, 0)
            masked = assoc_np * (adj_np != 0)

            # Build graph
            G = nx.from_numpy_array(masked)
            labels = {i: n for i, n in enumerate(assoc.columns)}
            nx.relabel_nodes(G, labels, copy=False)
            G.remove_edges_from(
                [(u, v) for u, v, d in G.edges(data=True) if d.get("weight", 0) == 0]
            )

            n_nodes = G.number_of_nodes()
            n_edges = G.number_of_edges()

            # Derive degree for node sizing; fall back to uniform size if unavailable
            cent = entry["centralities"]
            if cent is not None:
                cent_idx = cent.set_index("taxon")
                degrees = cent_idx["degree"].reindex(list(G.nodes()), fill_value=0)
                max_deg = degrees.max()
                node_sizes = (
                    (degrees / max_deg * 180 + 20).clip(lower=10).values
                    if max_deg > 0
                    else np.full(len(G.nodes()), 30)
                )
                top_nodes = set(cent_idx["degree"].nlargest(15).index) & set(G.nodes())
            else:
                node_sizes = np.full(len(G.nodes()), 30)
                top_nodes = set()
            node_labels = {n: n for n in G.nodes() if n in top_nodes}

            # Edge attributes
            weights = np.array([d["weight"] for _, _, d in G.edges(data=True)])
            edge_colors = [POS_COLOR if w > 0 else NEG_COLOR for w in weights]
            edge_widths = np.clip(np.abs(weights) * 2, 0.3, 3.0)

            # Spring layout (handles negative edge weights; Kamada-Kawai does not)
            np.random.seed(42)
            pos = nx.spring_layout(G, seed=42, k=1.5 / np.sqrt(max(n_nodes, 1)))

            fig, ax = plt.subplots(figsize=(16, 13))
            nx.draw_networkx_nodes(
                G,
                pos,
                node_size=node_sizes,
                node_color="#90CAF9",
                alpha=0.85,
                linewidths=0.4,
                edgecolors="white",
                ax=ax,
            )
            if n_edges > 0:
                nx.draw_networkx_edges(
                    G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.55, ax=ax
                )
            nx.draw_networkx_labels(
                G, pos, labels=node_labels, font_size=6, font_color="black", ax=ax
            )

            # Components: outline LCC nodes
            components = sorted(nx.connected_components(G), key=len, reverse=True)
            if len(components) > 1:
                lcc_nodes = list(components[0])
                lcc_pos = {n: pos[n] for n in lcc_nodes}
                lcc_sizes = (
                    node_sizes[[list(G.nodes()).index(n) for n in lcc_nodes]]
                    if len(G.nodes()) > 0
                    else []
                )
                nx.draw_networkx_nodes(
                    G.subgraph(lcc_nodes),
                    lcc_pos,
                    node_size=lcc_sizes,
                    node_color="#90CAF9",
                    alpha=0.0,
                    linewidths=1.5,
                    edgecolors="#1565C0",
                    ax=ax,
                )

            ax.set_title(
                f"{mth.upper()} — {lvl} level\n"
                f"{n_nodes} taxa, {n_edges} edges, "
                f"{len(components)} component(s)  |  "
                f"Node size ∝ degree  |  Labels = top-{min(15, len(top_nodes))} hubs",
                fontsize=11,
                fontweight="bold",
            )
            ax.axis("off")

            leg = [
                mpatches.Patch(
                    color=POS_COLOR, alpha=0.7, label="Positive association"
                ),
                mpatches.Patch(
                    color=NEG_COLOR, alpha=0.7, label="Negative association"
                ),
                mpatches.Patch(color="#90CAF9", label="Taxon node (size ∝ degree)"),
            ]
            ax.legend(
                handles=leg,
                loc="lower left",
                fontsize=9,
                framealpha=0.85,
                edgecolor="grey",
            )

            fig.tight_layout()
            save_page(pdf, fig)

    # ── Fig 5: Connected component size distributions ────────────────────
    print("  Fig 5: connected components")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        fig, axes = plt.subplots(
            1, len(methods_here), figsize=(5 * len(methods_here), 5), sharey=False
        )
        if len(methods_here) == 1:
            axes = [axes]
        fig.suptitle(
            f"Connected Component Size Distribution — {lvl}",
            fontsize=13,
            fontweight="bold",
        )

        for ax, mth in zip(axes, methods_here):
            entry = data[(lvl, mth)]
            assoc = entry["association_matrix"]
            adj = entry["adjacency_matrix"]
            if assoc is None or adj is None:
                ax.axis("off")
                continue

            assoc_np = assoc.values.copy()
            adj_np = adj.values.copy()
            np.fill_diagonal(assoc_np, 0)
            np.fill_diagonal(adj_np, 0)
            masked = assoc_np * (adj_np != 0)

            G = nx.from_numpy_array(masked)
            labels = {i: n for i, n in enumerate(assoc.columns)}
            nx.relabel_nodes(G, labels, copy=False)
            G.remove_edges_from(
                [(u, v) for u, v, d in G.edges(data=True) if d.get("weight", 0) == 0]
            )

            comp_sizes = sorted(
                [len(c) for c in nx.connected_components(G)], reverse=True
            )
            n_comps = len(comp_sizes)

            if n_comps <= 20:
                ax.bar(
                    range(1, n_comps + 1),
                    comp_sizes,
                    color=method_color(mth),
                    alpha=0.85,
                )
                ax.set_xlabel("Component rank")
            else:
                # Histogram for many small components
                ax.hist(
                    comp_sizes,
                    bins=min(30, len(set(comp_sizes))),
                    color=method_color(mth),
                    alpha=0.85,
                    edgecolor="white",
                )
                ax.set_xlabel("Component size (taxa)")

            ax.set_ylabel("Number of taxa")
            ax.set_title(
                f"{mth}\n{n_comps} component(s), "
                f"LCC={comp_sizes[0] if comp_sizes else 0} taxa",
                fontsize=9,
                fontweight="bold",
            )

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 6: Cluster / community composition ───────────────────────────
    print("  Fig 6: cluster compositions")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        # Determine grid layout
        n = len(methods_here)
        ncols = min(n, 3)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
        axes_flat = np.array(axes).flatten() if n > 1 else [axes]
        fig.suptitle(
            f"Cluster / Community Composition — {lvl}\n"
            "(clusters from fast_greedy algorithm on retained edges)",
            fontsize=13,
            fontweight="bold",
        )

        for ax, mth in zip(axes_flat, methods_here):
            clust = data[(lvl, mth)]["clusters"]
            if clust is None:
                ax.axis("off")
                continue

            # cluster 0 typically = unassigned/isolated; keep it for transparency
            cnt = clust["cluster"].value_counts().sort_index()
            n_taxa_total = len(clust)

            colors = plt.cm.Set2(np.linspace(0, 1, len(cnt)))
            bars = ax.bar(
                [f"C{c}" for c in cnt.index],
                cnt.values,
                color=colors,
                alpha=0.88,
                edgecolor="white",
            )
            for bar, val in zip(bars, cnt.values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    str(val),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

            ax.set_xlabel("Cluster")
            ax.set_ylabel("Number of taxa")
            ax.set_title(
                f"{mth}\n{cnt.nunique()} clusters, {n_taxa_total} taxa",
                fontsize=9,
                fontweight="bold",
            )
            ax.tick_params(axis="x", labelsize=8)

        # Hide unused subplots
        for ax in axes_flat[len(methods_here) :]:
            ax.axis("off")

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 7: Top hub taxa by degree ────────────────────────────────────
    print("  Fig 7: top hub taxa by degree")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        n = len(methods_here)
        fig, axes = plt.subplots(n, 1, figsize=(14, 4 * n))
        if n == 1:
            axes = [axes]
        fig.suptitle(
            f"Top {TOP_N} Taxa by Degree Centrality — {lvl}",
            fontsize=13,
            fontweight="bold",
        )

        for ax, mth in zip(axes, methods_here):
            cent = data[(lvl, mth)]["centralities"]
            clust = data[(lvl, mth)]["clusters"]
            if cent is None:
                ax.axis("off")
                continue

            top = cent.nlargest(TOP_N, "degree").copy()
            if clust is not None:
                top = top.merge(clust[["taxon", "cluster"]], on="taxon", how="left")
                top["cluster"] = top["cluster"].fillna(-1).astype(int)
                n_clusters = top["cluster"].nunique()
                cmap = plt.cm.Set2(np.linspace(0, 1, max(n_clusters, 2)))
                cluster_ids = sorted(top["cluster"].unique())
                cmap_dict = {c: cmap[i % len(cmap)] for i, c in enumerate(cluster_ids)}
                bar_colors = [cmap_dict[c] for c in top["cluster"]]
            else:
                bar_colors = method_color(mth)

            ax.barh(
                range(len(top)),
                top["degree"][::-1].values,
                color=bar_colors[::-1] if isinstance(bar_colors, list) else bar_colors,
                alpha=0.88,
            )
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(top["taxon"].iloc[::-1].values, fontsize=7)
            ax.set_xlabel("Degree centrality")
            ax.set_title(f"{mth}", fontsize=10, fontweight="bold")

            if clust is not None and top["cluster"].nunique() > 1:
                legend_handles = [
                    mpatches.Patch(color=cmap_dict[c], label=f"Cluster {c}")
                    for c in cluster_ids
                ]
                ax.legend(
                    handles=legend_handles,
                    fontsize=7,
                    loc="lower right",
                    ncol=min(6, len(cluster_ids)),
                )

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 8: Top hub taxa by betweenness ───────────────────────────────
    print("  Fig 8: top hub taxa by betweenness")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        n = len(methods_here)
        fig, axes = plt.subplots(n, 1, figsize=(14, 4 * n))
        if n == 1:
            axes = [axes]
        fig.suptitle(
            f"Top {TOP_N} Taxa by Betweenness Centrality — {lvl}\n"
            "(high betweenness = bridge between communities)",
            fontsize=13,
            fontweight="bold",
        )

        for ax, mth in zip(axes, methods_here):
            cent = data[(lvl, mth)]["centralities"]
            bmax = cent["betweenness"].max() if cent is not None else None
            if cent is None or bmax != bmax or bmax == 0:  # None, NaN, or zero
                ax.axis("off")
                ax.set_title(f"{mth} — no betweenness data", fontsize=9)
                continue

            top = cent[cent["betweenness"] > 0].nlargest(TOP_N, "betweenness").copy()
            ax.barh(
                range(len(top)),
                top["betweenness"][::-1].values,
                color=method_color(mth),
                alpha=0.88,
            )
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(top["taxon"].iloc[::-1].values, fontsize=7)
            ax.set_xlabel("Betweenness centrality")
            ax.set_title(
                f"{mth}  ({len(top)} taxa with betweenness > 0)",
                fontsize=10,
                fontweight="bold",
            )

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 9: Degree distributions ──────────────────────────────────────
    print("  Fig 9: degree distributions")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        fig, axes = plt.subplots(
            1, len(methods_here), figsize=(5 * len(methods_here), 5)
        )
        if len(methods_here) == 1:
            axes = [axes]
        fig.suptitle(f"Degree Distribution — {lvl}", fontsize=13, fontweight="bold")

        for ax, mth in zip(axes, methods_here):
            cent = data[(lvl, mth)]["centralities"]
            if cent is None:
                ax.axis("off")
                continue
            deg = cent["degree"]
            ax.hist(
                deg, bins=30, color=method_color(mth), alpha=0.85, edgecolor="white"
            )
            ax.axvline(
                deg.median(),
                color="black",
                linestyle="--",
                linewidth=1.2,
                label=f"Median={deg.median():.2f}",
            )
            ax.set_xlabel("Degree centrality")
            ax.set_ylabel("Number of taxa")
            ax.set_title(mth, fontsize=10, fontweight="bold")
            ax.legend(fontsize=8)

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 10: Top taxa by closeness centrality ─────────────────────────
    print("  Fig 10: closeness centrality")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        n = len(methods_here)
        fig, axes = plt.subplots(n, 1, figsize=(14, 4 * n))
        if n == 1:
            axes = [axes]
        fig.suptitle(
            f"Top {TOP_N} Taxa by Closeness Centrality — {lvl}\n"
            "(high closeness = shortest average distance to all other taxa)",
            fontsize=13,
            fontweight="bold",
        )

        for ax, mth in zip(axes, methods_here):
            cent = data[(lvl, mth)]["centralities"]
            cmax = cent["closeness"].max() if cent is not None else None
            if cent is None or cmax != cmax or cmax == 0:  # None, NaN, or zero
                ax.axis("off")
                ax.set_title(f"{mth} — no closeness data", fontsize=9)
                continue

            top = cent[cent["closeness"] > 0].nlargest(TOP_N, "closeness").copy()
            ax.barh(
                range(len(top)),
                top["closeness"][::-1].values,
                color=method_color(mth),
                alpha=0.88,
            )
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(top["taxon"].iloc[::-1].values, fontsize=7)
            ax.set_xlabel("Closeness centrality")
            ax.set_title(mth, fontsize=10, fontweight="bold")

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 11: Positive vs negative edge breakdown ───────────────────────
    print("  Fig 11: positive vs negative edges")

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if not methods_here:
            continue

        rows = []
        for mth in methods_here:
            assoc = data[(lvl, mth)]["association_matrix"]
            adj = data[(lvl, mth)]["adjacency_matrix"]
            if assoc is None or adj is None:
                continue
            a = assoc.values.copy()
            b = adj.values.copy()
            np.fill_diagonal(a, 0)
            np.fill_diagonal(b, 0)
            masked = a * (b != 0)
            # count upper triangle only
            ut = np.triu(masked, k=1)
            n_pos = int((ut > 0).sum())
            n_neg = int((ut < 0).sum())
            rows.append({"method": mth, "Positive": n_pos, "Negative": abs(n_neg)})

        if not rows:
            continue

        pe_df = pd.DataFrame(rows).set_index("method")

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(
            f"Positive vs Negative Edges — {lvl}", fontsize=13, fontweight="bold"
        )

        # Stacked bar
        ax = axes[0]
        x = np.arange(len(pe_df))
        ax.bar(x, pe_df["Positive"], color=POS_COLOR, alpha=0.85, label="Positive")
        ax.bar(
            x,
            pe_df["Negative"],
            bottom=pe_df["Positive"],
            color=NEG_COLOR,
            alpha=0.85,
            label="Negative",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(pe_df.index, rotation=25, ha="right")
        ax.set_ylabel("Edge count")
        ax.set_title("Raw edge counts", fontsize=10, fontweight="bold")
        ax.legend(fontsize=9)

        # Proportion bar
        ax = axes[1]
        totals = pe_df["Positive"] + pe_df["Negative"]
        pct_pos = pe_df["Positive"] / totals * 100
        pct_neg = pe_df["Negative"] / totals * 100
        ax.bar(x, pct_pos, color=POS_COLOR, alpha=0.85, label="Positive")
        ax.bar(
            x, pct_neg, bottom=pct_pos, color=NEG_COLOR, alpha=0.85, label="Negative"
        )
        ax.axhline(50, color="grey", linestyle="--", linewidth=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(pe_df.index, rotation=25, ha="right")
        ax.set_ylabel("% of edges")
        ax.set_ylim(0, 110)
        ax.set_title("Edge proportion", fontsize=10, fontweight="bold")
        ax.legend(fontsize=9)

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 12 & 13: Method agreement heatmaps ───────────────────────────
    print("  Fig 12-13: method agreement heatmaps")

    def plot_heatmap_grid(
        comp_dict_key, title, subtitle, cmap, vmin, vmax, center=None
    ):
        mats = {
            lvl: comp_data[lvl][comp_dict_key]
            for lvl in TAX_LEVELS
            if comp_data[lvl][comp_dict_key] is not None
        }
        if not mats:
            return
        n = len(mats)
        ncols = min(n, 2)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows))
        axes_flat = np.array(axes).flatten() if n > 1 else [axes]
        fig.suptitle(f"{title}\n{subtitle}", fontsize=13, fontweight="bold")
        for ax, (lvl, mat) in zip(axes_flat, mats.items()):
            mask = np.triu(np.ones_like(mat.values, dtype=bool), k=1)
            sns.heatmap(
                mat,
                mask=mask,
                annot=True,
                fmt=".2f",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                center=center,
                linewidths=0.5,
                linecolor="white",
                ax=ax,
                annot_kws={"size": 8},
            )
            ax.set_title(lvl, fontsize=11, fontweight="bold")
            ax.tick_params(axis="x", rotation=35, labelsize=8)
            ax.tick_params(axis="y", rotation=0, labelsize=8)
        for ax in axes_flat[len(mats) :]:
            ax.axis("off")
        fig.tight_layout()
        save_page(pdf, fig)

    plot_heatmap_grid(
        "edge_overlap",
        "Method Agreement: Edge Jaccard Overlap",
        "1.0 = identical edge sets   0.0 = no shared edges",
        cmap="Blues",
        vmin=0,
        vmax=1,
    )
    plot_heatmap_grid(
        "deg_corr",
        "Method Agreement: Degree Centrality Correlation (Spearman r)",
        "Positive = methods rank taxa similarly   Negative = opposite rankings",
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        center=0,
    )
    plot_heatmap_grid(
        "hub_overlap",
        "Hub Taxa Agreement Across Methods (Jaccard)",
        "Which methods identify the same taxa as network hubs?",
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
    )

    # ── Fig 14: Shared hub taxa across methods (upset-style) ─────────────
    print("  Fig 14: shared hub taxa across methods")

    HUB_TOP_N = 20  # define hubs as top-N by degree within each method

    for lvl in TAX_LEVELS:
        methods_here = [m for m in present_methods if (lvl, m) in data]
        if len(methods_here) < 2:
            continue

        hub_sets = {}
        for mth in methods_here:
            cent = data[(lvl, mth)]["centralities"]
            if cent is None:
                continue
            top = cent.nlargest(HUB_TOP_N, "degree")["taxon"].tolist()
            hub_sets[mth] = set(top)

        if len(hub_sets) < 2:
            continue

        # All taxa appearing in any hub set
        all_hubs = sorted(set.union(*hub_sets.values()))

        # Build presence matrix
        pres = pd.DataFrame(
            {mth: [t in hub_sets[mth] for t in all_hubs] for mth in hub_sets},
            index=all_hubs,
        )
        pres["n_methods"] = pres.sum(axis=1)
        pres = pres.sort_values("n_methods", ascending=False)

        fig, axes = plt.subplots(
            1,
            2,
            figsize=(16, max(6, len(all_hubs) * 0.28 + 2)),
            gridspec_kw={"width_ratios": [3, 1]},
        )
        fig.suptitle(
            f"Hub Taxa Agreement — {lvl}\n(top-{HUB_TOP_N} taxa by degree per method)",
            fontsize=13,
            fontweight="bold",
        )

        # Heatmap of presence/absence
        ax = axes[0]
        sns.heatmap(
            pres[list(hub_sets.keys())].astype(int),
            annot=False,
            cmap="Blues",
            vmin=0,
            vmax=1,
            linewidths=0.3,
            linecolor="white",
            ax=ax,
            cbar=False,
        )
        ax.set_xlabel("Method")
        ax.set_ylabel("Taxon")
        ax.tick_params(axis="x", rotation=35, labelsize=9)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_title("Presence in top-hub set", fontsize=10, fontweight="bold")

        # Bar: n_methods
        ax2 = axes[1]
        colors_bar = plt.cm.YlOrRd(pres["n_methods"] / len(hub_sets))
        ax2.barh(range(len(pres)), pres["n_methods"], color=colors_bar, alpha=0.9)
        ax2.set_yticks(range(len(pres)))
        ax2.set_yticklabels(pres.index, fontsize=7)
        ax2.set_xlabel("Methods agreeing")
        ax2.set_title("Agreement count", fontsize=10, fontweight="bold")
        ax2.set_xlim(0, len(hub_sets) + 0.5)
        ax2.axvline(len(hub_sets), color="grey", linestyle="--", linewidth=0.8)

        fig.tight_layout()
        save_page(pdf, fig)

    # ── Fig 15: Key findings summary ─────────────────────────────────────
    print("  Fig 15: key findings summary")

    # Compute a few auto-populated stats for the summary
    summary_lines = [
        "KEY FINDINGS SUMMARY",
        "─" * 65,
        "",
        "DATA",
        f"  Taxonomic level(s): {', '.join(TAX_LEVELS)}",
        f"  Methods run:        {', '.join(present_methods)}",
        "",
        "NETWORK TOPOLOGY",
    ]

    for lvl in TAX_LEVELS:
        cd = comp_data[lvl]["basic_stats"]
        if cd is None:
            continue
        summary_lines.append(f"\n  [{lvl}]")
        for _, row in cd.iterrows():
            mth = row["method"]
            if mth not in present_methods:
                continue
            density = row.get("edge_density", float("nan"))
            lcc = row.get("lcc_size", float("nan"))
            ncomp = row.get("n_components", float("nan"))
            n_edges = row.get("n_edges", float("nan"))
            lcc_str = f"{int(lcc):>4}" if lcc == lcc else "  NA"
            ncomp_str = str(int(ncomp)) if ncomp == ncomp else "NA"
            edges_str = f"{int(n_edges):>6}" if n_edges == n_edges else "    NA"
            dens_str = f"{density:.3f}" if density == density else "   NA"
            summary_lines.append(
                f"    {mth:12s}  edges={edges_str}  "
                f"density={dens_str}  LCC={lcc_str} taxa  "
                f"components={ncomp_str}"
            )

    summary_lines += [
        "",
        "INTERPRETATION NOTES",
        "  • Sparse methods (SPRING, SpiecEasi, CCLasso) retain only statistically",
        "    supported edges; dense methods (SparCC, Spearman) do not sparsify.",
        "  • Nodes = taxa; edges = significant co-occurrence associations.",
        "  • Positive edge: taxa tend to co-occur across samples.",
        "  • Negative edge: taxa are mutually exclusive across samples.",
        "  • Hub taxa (high degree/betweenness) are potential keystone taxa.",
        "  • Connected components = groups of taxa with no cross-group associations.",
        "  • Cluster 0 (fast_greedy) = unassigned/isolated taxa.",
        "",
        "CAUTION",
        "  • Network structure depends strongly on the association method.",
        "  • Compare hub taxa across methods to identify robust findings.",
        "  • Small sample sizes reduce statistical power for edge detection.",
    ]

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.axis("off")
    ax.text(
        0.02,
        0.98,
        "\n".join(summary_lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#F5F5F5", alpha=0.8),
    )
    fig.tight_layout()
    save_page(pdf, fig)

print(f"\nDone! Output: {OUTPUT_PDF}")
print("=" * 55)
