from __future__ import annotations
import json
from typing import Any
import zarr
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from cherita.utils.adata_utils import (
    continuous2categorical,
    get_group_index,
    get_indices_in_array,
    parse_data,
)

MARKER_PIXEL_SIZE = 50


def dotplot(
    adata_group: zarr.Group,
    markers: list[str],
    obs_col: dict,
    expression_cutoff: float = 0.0,
    mean_only_expressed: bool = False,
    standard_scale: str = None,
) -> Any:
    """Method to generate a Plotly dotplot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        markers (list[str]): Root zarr Group of an Anndata-Zarr object
        obs_col (str): The obs column to group data
        expression_cutoff (float, optional): Minimum expression value to consider
            if mean_only_expressed is set to True. Defaults to 0.0.
        mean_only_expressed (bool, optional): Whether to only average values
            above the expression cutoff. Defaults to False.
        standard_scale (str, optional): Scaling method to use.
            Can be set to None, "var" or "group. Defaults to None.

    Returns:
        Any: A Plotly dotplot JSON as a Python object
    """
    marker_idx = get_indices_in_array(get_group_index(adata_group.var), markers)
    obs_colname = obs_col["name"]
    obs = parse_data(adata_group.obs[obs_colname])

    df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=markers)

    if obs_col["type"] == "continuous":
        df[obs_colname] = continuous2categorical(obs, obs_col["bins"]["thresholds"])

    elif obs_col["type"] == "categorical":
        df[obs_colname] = obs

    df.set_index(obs_colname, append=True, inplace=True)

    # Following Scanpy's dotplot
    # https://github.com/scverse/scanpy/blob/ed3b277b2f498e3cab04c9416aaddf97eec8c3e2/scanpy/plotting/_dotplot.py

    bool_df = df > expression_cutoff

    dot_size_df = (
        bool_df.groupby(obs_colname).sum() / bool_df.groupby(obs_colname).count()
    )

    if mean_only_expressed:
        dot_color_df = df.mask(~bool_df).groupby(obs_colname).mean().fillna(0)
    else:
        dot_color_df = df.groupby(obs_colname).mean()
    mean_df = dot_color_df

    if standard_scale == "group":
        dot_color_df = dot_color_df.sub(dot_color_df.min(1), axis=0)
        dot_color_df = dot_color_df.div(dot_color_df.max(1), axis=0).fillna(0)
    elif standard_scale == "var":
        dot_color_df -= dot_color_df.min(0)
        dot_color_df = (dot_color_df / dot_color_df.max(0)).fillna(0)
    elif standard_scale is None:
        pass

    x = df.index.get_level_values(obs_colname).categories

    data = [
        go.Scatter(
            x=x,
            y=[col] * len(x),
            name=col,
            mode="markers",
            marker=dict(
                color=dot_color_df[col],
                coloraxis="coloraxis",
                line=dict(width=0.5, color="black"),
                size=dot_size_df[col] * MARKER_PIXEL_SIZE,
            ),
            customdata=np.dstack((mean_df[col], dot_size_df[col]))[0],
            hovertemplate="Mean expression: %{customdata[0]:.3f} <br>"
            "Fraction of cells in group: %{customdata[1]:.3f}",
            showlegend=False,
        )
        for col in df.columns
    ]

    marker_legends = np.arange(0.2, 1.2, 0.2)
    marker_legends_data = go.Scatter(
        x=[1] * len(marker_legends),
        y=marker_legends * 100,
        mode="markers",
        marker=dict(
            color="gray",
            size=marker_legends * MARKER_PIXEL_SIZE,
            line=dict(width=0.5, color="black"),
        ),
        showlegend=False,
        hovertemplate=None,
        hoverinfo="none",
    )

    fig = make_subplots(
        rows=1, cols=2, column_widths=[0.9, 0.1], horizontal_spacing=0.10
    )

    for d in data:
        fig.add_trace(d, row=1, col=1)
    fig.add_trace(marker_legends_data, row=1, col=2)
    fig.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis_type="category",
        coloraxis=dict(
            colorscale=None,
            colorbar_x=0.82,
            colorbar=dict(title=dict(text="Mean expression in group", side="right")),
        ),
        xaxis=dict(
            title=obs_colname
            + (
                " ({} bins)".format(obs_col["bins"]["nBins"])
                if obs_col["type"] == "continuous"
                else ""
            ),
            showline=True,
            linewidth=1,
            linecolor="black",
            mirror=True,
        ),
        yaxis=dict(
            title="Markers", showline=True, linewidth=1, linecolor="black", mirror=True
        ),
        xaxis2=dict(
            zeroline=False,
            showline=False,
            showticklabels=False,
            fixedrange=True,
            showgrid=False,
        ),
        yaxis2=dict(
            title="Fraction of cells in group (%)",
            zeroline=False,
            showline=False,
            showgrid=False,
            fixedrange=True,
            side="right",
            tickvals=marker_legends * 100,
        ),
    )

    data_values_range = {
        "min": float(dot_color_df.values.min()),
        "max": float(dot_color_df.values.max()),
    }

    response_obj = json.loads(fig.to_json())
    response_obj["range"] = data_values_range

    return response_obj
