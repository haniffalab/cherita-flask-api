from __future__ import annotations
import json
from typing import Any, Union
import zarr
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from cherita.resources.errors import BadRequest, InvalidObs

from cherita.utils.adata_utils import (
    to_categorical,
    parse_data,
)
from cherita.utils.models import Marker


def matrixplot(
    adata_group: zarr.Group,
    var_keys: list[Union[int, str, dict]],
    obs_col: dict,
    obs_values: list[str] = None,
    standard_scale: str = None,
    var_names_col: str = None,
    obs_indices: list[int] = None,
) -> Any:
    """Method to generate a Plotly matrixplot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        var_keys (list[Union[int, str, dict]]): List of markers present in var.
        obs_col (dict): The obs column to group data.
        obs_values (list[str], optional): List of values in obs to plot.
        standard_scale (str, optional): Scaling method to use.
            Can be set to None, "var" or "group. Defaults to None.
        var_names_col (str, optional): Column in var to pull markers' names from.
            Defaults to None.

    Returns:
        Any: A Plotly matrixplot JSON as a Python object
    """
    if not isinstance(obs_col, dict):
        raise BadRequest("'selectedObs' must be an object")

    markers = [Marker.from_any(adata_group, v, var_names_col) for v in var_keys]

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

    if obs_indices is not None:
        df = df.iloc[obs_indices]

    if obs_values is not None:
        df = df[df[obs_colname].isin(obs_values)]
        df[obs_colname] = df[obs_colname].cat.remove_unused_categories()

    values_df = df.groupby(obs_colname, observed=False).mean()

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
            coloraxis="coloraxis",
        ),
        layout=dict(
            autosize=False,
            plot_bgcolor="rgba(0,0,0,0)",
            yaxis_type="category",
            xaxis_type="category",
            xaxis=dict(title="Markers", tickvals=values_df.columns),
            yaxis=dict(
                title=obs_colname + (f" ({bins} bins)" if bins else ""),
                tickvals=values_df.index,
                scaleanchor="x",
            ),
        ),
    )

    return json.loads(fig.to_json())
