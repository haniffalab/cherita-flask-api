from __future__ import annotations
import math
import json
from typing import Any
import zarr
import pandas as pd
import plotly.graph_objects as go
from cherita.resources.errors import BadRequest, InvalidKey, InvalidObs, InvalidVar

from cherita.utils.adata_utils import (
    get_group_index,
    get_indices_in_array,
    parse_data,
    to_categorical,
)

CHUNK_SIZE = 60000


def split_df(df, chunk_size):
    chunks = []
    n_chunks = math.ceil(len(df) / chunk_size)
    for i in range(n_chunks):
        chunks.append(
            # To address gaps between image html elements
            # Include previous items to force overlap
            df[
                i * chunk_size
                - (int(len(df) * 0.001) if i > 0 else 0) : (i + 1) * chunk_size
            ]
        )
    return chunks


def heatmap(adata_group: zarr.Group, markers: list[str], obs_col: dict) -> Any:
    """Method to generate a Plotly heatmap plot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        markers (list[str]): List of markers present in var.
        obs_col (dict): The obs column to group data

    Returns:
        Any: A Plotly heatmap plot JSON as a Python object
    """
    if not isinstance(obs_col, dict):
        raise BadRequest("'selectedObs' must be an object")

    try:
        marker_idx = get_indices_in_array(get_group_index(adata_group.var), markers)
    except InvalidKey:
        raise InvalidVar(f"Invalid features {markers}")

    obs_colname = obs_col["name"]
    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=markers)

    layout = dict(yaxis=dict(title="Markers"))

    df[obs_colname], bins = to_categorical(obs, **obs_col)

    df = df.sort_values(by=[obs_colname])
    df = df.reset_index()

    ticks = set([])
    middle_ticks = []
    for c in df[obs_colname].cat.categories:
        ticks.add(df[df[obs_colname] == c][markers].index.min())
        ticks.add(df[df[obs_colname] == c][markers].index.max() + 1)
        middle_ticks.append(
            (
                df[df[obs_colname] == c][markers].index.min()
                + (
                    df[df[obs_colname] == c][markers].index.max()
                    - df[df[obs_colname] == c][markers].index.min()
                )
                // 2
            )
        )
    ticks = list(ticks)
    ticks.sort()

    layout.update(
        dict(
            xaxis=dict(
                title=obs_colname + (f" ({bins} bins)" if bins else ""),
                tickvals=middle_ticks,
                ticktext=list(df[obs_colname].cat.categories),
                minor=dict(tickvals=ticks, ticks="outside", ticklen=5),
            )
        )
    )

    # To handle data above 65k rows (image limit)
    sub_dfs = split_df(df, CHUNK_SIZE)
    sub_heatmaps = []
    for d in sub_dfs:
        sub_heatmaps.append(
            go.Heatmap(
                z=d[markers].transpose(),
                y=markers,
                x=d.index.union([d.index.max() + 1]),
                coloraxis="coloraxis",
                name="",
            )
        )

    fig = go.Figure(data=sub_heatmaps, layout=layout)

    if not len(markers) > 2:
        fig.update_yaxes(fixedrange=True)

    return json.loads(fig.to_json())
