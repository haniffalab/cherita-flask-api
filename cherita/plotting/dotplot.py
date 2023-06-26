import json
import zarr
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.validator_cache import ValidatorCache

from cherita.utils.adata_utils import get_group_index, get_indices_in_array, parse_data

MARKER_PIXEL_SIZE = 50


def dotplot(
    adata_group: zarr.Group,
    markers: list[str] = None,
    obs_col: str = None,
    expression_cutoff=0.5,
    mean_only_expressed=False,
    standard_scale=None,
):
    marker_idx = get_indices_in_array(get_group_index(adata_group.var), markers)
    obs = parse_data(adata_group.obs[obs_col])

    df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=markers)

    df[obs_col] = obs
    df.set_index(obs_col, append=True, inplace=True)

    # Following Scanpy's dotplot
    # https://github.com/scverse/scanpy/blob/ed3b277b2f498e3cab04c9416aaddf97eec8c3e2/scanpy/plotting/_dotplot.py

    bool_df = df > expression_cutoff

    dot_size_df = bool_df.groupby(obs_col).sum() / bool_df.groupby(obs_col).count()

    if mean_only_expressed:
        dot_color_df = df.mask(~bool_df).groupby(obs_col).mean().fillna(0)
    else:
        dot_color_df = df.groupby(obs_col).mean()
    mean_df = dot_color_df

    if standard_scale == "group":
        dot_color_df = dot_color_df.sub(dot_color_df.min(1), axis=0)
        dot_color_df = dot_color_df.div(dot_color_df.max(1), axis=0).fillna(0)
    elif standard_scale == "var":
        dot_color_df -= dot_color_df.min(0)
        dot_color_df = (dot_color_df / dot_color_df.max(0)).fillna(0)
    elif standard_scale is None:
        pass

    x = df.index.get_level_values(obs_col).categories

    data = [
        go.Scatter(
            x=x,
            y=[col] * len(x),
            mode="markers",
            marker=dict(
                color=dot_color_df[col],
                coloraxis="coloraxis",
                line=dict(width=0.5, color="black"),
                size=dot_size_df[col] * MARKER_PIXEL_SIZE,
            ),
            customdata=np.dstack((mean_df[col], dot_size_df[col]))[0],
            hovertemplate="Mean expression: %{customdata[0]:.3f} <br> Fraction of cells in group: %{customdata[1]:.3f}",
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
        coloraxis=dict(colorscale=None, colorbar_x=0.82),
        xaxis=dict(
            title=obs_col, showline=True, linewidth=1, linecolor="black", mirror=True
        ),
        yaxis=dict(
            title="Markers", showline=True, linewidth=1, linecolor="black", mirror=True
        ),
        xaxis2=dict(
            title="Fraction of cells in group (%)",
            zeroline=False,
            showline=False,
            showticklabels=False,
            fixedrange=True,
            showgrid=False,
            title_standoff=25,
        ),
        yaxis2=dict(
            zeroline=False,
            showline=False,
            showgrid=False,
            fixedrange=True,
            side="right",
            tickvals=marker_legends * 100,
        ),
    )

    return json.loads(fig.to_json())
