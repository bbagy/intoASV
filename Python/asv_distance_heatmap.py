#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
asv_distance_heatmap.py (FINAL — no ASV in hover)

Usage:
  python3 asv_distance_heatmap.py \
      --dist  asv_distance_matrix_251204.csv \
      --anno  cluster_diagnostics_251204.csv \
      --order asv_order_251204.csv \
      --out   ASV_distance_heatmap_251204.html
"""

import argparse
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create interactive ASV distance heatmap (ClusterID + Species only)."
    )
    parser.add_argument("--dist", required=True,
                        help="ASV × ASV distance matrix CSV")
    parser.add_argument("--anno", required=True,
                        help="cluster_diagnostics CSV (ASV, ClusterID, Species)")
    parser.add_argument("--order", required=True,
                        help="CSV file with ASV order (one ASV per line)")
    parser.add_argument("--out", required=True,
                        help="Output HTML")

    return parser.parse_args()


def main():
    args = parse_args()

    # -----------------------------
    # 1) Load distance matrix
    # -----------------------------
    dist_df = pd.read_csv(args.dist, index_col=0)
    asvs_matrix = list(dist_df.index)

    # -----------------------------
    # 2) Load ASV order from R
    # -----------------------------
    order_df = pd.read_csv(args.order, header=None)
    asv_order = order_df[0].astype(str).tolist()

    # keep only ASVs present in matrix
    asv_order = [a for a in asv_order if a in asvs_matrix]

    # enforce order
    dist_df = dist_df.loc[asv_order, asv_order]
    dist = dist_df.values

    # -----------------------------
    # 3) Load annotation (cluster diagnostics)
    # -----------------------------
    anno = pd.read_csv(args.anno)

    required = {"ASV", "ClusterID", "Species"}
    if not required.issubset(anno.columns):
        missing = required - set(anno.columns)
        raise ValueError(f"Annotation missing: {missing}")

    # filter & reorder annotation
    anno = anno[anno["ASV"].astype(str).isin(asv_order)]
    anno = anno.drop_duplicates(subset="ASV")
    anno = anno.set_index("ASV").loc[asv_order]

    clusters = anno["ClusterID"].astype(str).to_numpy()
    species  = anno["Species"].astype(str).fillna("NA").to_numpy()

    # -----------------------------
    # 4) Build hover customdata
    # -----------------------------
    n = len(asv_order)

    row_cluster  = clusters.reshape(-1, 1)
    col_cluster  = clusters.reshape(1,  -1)
    row_species  = species.reshape(-1, 1)
    col_species  = species.reshape(1,  -1)

    custom = np.stack([
        np.repeat(row_cluster, n, axis=1),
        np.repeat(col_cluster, n, axis=0),
        np.repeat(row_species, n, axis=1),
        np.repeat(col_species, n, axis=0),
    ], axis=2)

    # -----------------------------
    # 5) Create heatmap
    # -----------------------------
    fig = go.Figure(
        go.Heatmap(
            z=dist,
            x=asv_order,
            y=asv_order,
            colorscale="YlOrRd",
            colorbar=dict(title="Distance"),
            customdata=custom,
            hovertemplate=
                "Row Cluster: %{customdata[0]}<br>"
                "Col Cluster: %{customdata[1]}<br>"
                "Row Species: %{customdata[2]}<br>"
                "Col Species: %{customdata[3]}<br>"
                "Distance: %{z:.4f}<extra></extra>",
        )
    )

    # hide unreadable axis tick labels
    fig.update_layout(
        title=dict(
            text="ASV-level Pairwise Distance Heatmap",
            x=0.5,
            xanchor='center',
            yanchor='top',
            font=dict(size=20)
        ),
        xaxis=dict(showticklabels=False),
        yaxis=dict(showticklabels=False),
        width=1200,
        height=1200
    )
    # -----------------------------
    # 6) Save HTML (PlotlyJS via CDN)
    # -----------------------------
    out_path = os.path.abspath(args.out)
    fig.write_html(out_path, include_plotlyjs="cdn")
    print(f"[asv_distance_heatmap] HTML saved → {out_path}")


if __name__ == "__main__":
    main()
