from __future__ import annotations
import json
import zarr
import pandas as pd
import plotly.graph_objects as go

from cherita.utils.adata_utils import get_group_index, get_indices_in_array, parse_data


def matrixplot(
    adata_group: zarr.Group,
    markers: list[str] = None,
    obs_col: str = None,
    standard_scale: str = None,
):
    marker_idx = get_indices_in_array(get_group_index(adata_group.var), markers)
    obs = parse_data(adata_group.obs[obs_col])

    df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=markers)
    df[obs_col] = obs

    values_df = df.groupby(obs_col).mean()

    if standard_scale == "group":
        values_df = values_df.sub(values_df.min(1), axis=0)
        values_df = values_df.div(values_df.max(1), axis=0).fillna(0)
    elif standard_scale == "var":
        values_df -= values_df.min(0)
        values_df = (values_df / values_df.max(0)).fillna(0)
    elif standard_scale is None:
        pass

    fig = go.Figure(
        data=go.Heatmap(
            z=values_df.values,
            x=values_df.columns,
            y=values_df.index,
            colorbar=dict(title=dict(text="Mean expression in group", side="right")),
        ),
        layout=dict(
            autosize=False,
            plot_bgcolor="rgba(0,0,0,0)",
            yaxis_type="category",
            xaxis_type="category",
            xaxis=dict(title="Markers"),
            yaxis=dict(title=obs_col, tickvals=values_df.index, scaleanchor="x"),
        ),
    )

    return json.loads(fig.to_json())
