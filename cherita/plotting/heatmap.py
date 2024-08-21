from __future__ import annotations
import math
import json
from typing import Union, Any
import zarr
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from cherita.resources.errors import BadRequest, InvalidObs

from cherita.utils.adata_utils import parse_data, to_categorical
from cherita.utils.models import Marker

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


def heatmap(
    adata_group: zarr.Group,
    var_keys: list[Union[int, str, dict]],
    obs_col: dict,
    obs_values: list[str] = None,
    var_names_col: str = None,
) -> Any:
    """Method to generate a Plotly heatmap plot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        var_keys (list[Union[int, str, dict]]): List of markers present in var.
        obs_col (dict): The obs column to group data
        obs_values (list[str], optional): List of values in obs to plot.
        var_names_col (str, optional): Column in var to pull markers' names from.
            Defaults to None.

    Returns:
        Any: A Plotly heatmap plot JSON as a Python object
    """
    if not isinstance(obs_col, dict):
        raise BadRequest("'selectedObs' must be an object")

    markers = [Marker.from_any(adata_group, v, var_names_col) for v in var_keys]
    marker_names = [m.name for m in markers]

    obs_colname = obs_col["name"]
    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    df = pd.DataFrame(
        np.concatenate([[m.X] for m in markers]).T,
        columns=[m.name for m in markers],
    )

    df[obs_colname], bins = to_categorical(obs, **obs_col)

    if obs_values is not None:
        df = df[df[obs_colname].isin(obs_values)]
        df[obs_colname] = df[obs_colname].cat.remove_unused_categories()

    df = df.sort_values(by=[obs_colname])
    df = df.reset_index()

    layout = dict(yaxis=dict(title="Markers"))
    ticks = set([])
    middle_ticks = []
    for c in df[obs_colname].cat.categories:
        ticks.add(df[df[obs_colname] == c][marker_names].index.min())
        ticks.add(df[df[obs_colname] == c][marker_names].index.max() + 1)
        middle_ticks.append(
            (
                df[df[obs_colname] == c][marker_names].index.min()
                + (
                    df[df[obs_colname] == c][marker_names].index.max()
                    - df[df[obs_colname] == c][marker_names].index.min()
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
                z=d[marker_names].transpose(),
                y=marker_names,
                x=d.index.union([d.index.max() + 1]),
                coloraxis="coloraxis",
                name="",
            )
        )

    fig = go.Figure(data=sub_heatmaps, layout=layout)

    if not len(var_keys) > 2:
        fig.update_yaxes(fixedrange=True)

    return json.loads(fig.to_json())
