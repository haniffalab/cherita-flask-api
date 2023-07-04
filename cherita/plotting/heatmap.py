from __future__ import annotations
import json
from typing import Any
import zarr
import pandas as pd
import plotly.graph_objects as go

from cherita.utils.adata_utils import get_group_index, get_indices_in_array, parse_data


def heatmap(adata_group: zarr.Group, markers: list[str], obs_col: str) -> Any:
    """Method to generate a Plotly heatmap plot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        markers (list[str]): List of markers present in var.
        obs_col (str): The obs column to use for the X axis

    Returns:
        Any: A Plotly heatmap plot JSON as a Python object
    """
    marker_idx = get_indices_in_array(get_group_index(adata_group.var), markers)
    obs = parse_data(adata_group.obs[obs_col])

    df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=markers)

    layout = dict(yaxis=dict(title="Markers"))

    if obs.dtype == "category":
        df[obs_col] = obs

        df = df.sort_values(by=[obs_col])
        df = df.reset_index()

        ticks = set([])
        middle_ticks = []
        for c in df[obs_col].cat.categories:
            ticks.add(df[df[obs_col] == c][markers].index.min())
            ticks.add(df[df[obs_col] == c][markers].index.max() + 1)
            middle_ticks.append(
                (
                    df[df[obs_col] == c][markers].index.min()
                    + (
                        df[df[obs_col] == c][markers].index.max()
                        - df[df[obs_col] == c][markers].index.min()
                    )
                    // 2
                )
            )
        ticks = list(ticks)
        ticks.sort()
        ticks[-1] -= 1

        layout.update(
            dict(
                xaxis=dict(
                    title=obs_col,
                    tickvals=middle_ticks,
                    ticktext=list(df[obs_col].cat.categories),
                    minor=dict(tickvals=ticks, ticks="outside", ticklen=5),
                )
            )
        )

    fig = go.Figure(go.Heatmap(z=df[markers].transpose(), y=markers), layout=layout)

    if not len(markers) > 2:
        fig.update_yaxes(fixedrange=True)

    return json.loads(fig.to_json())
