from __future__ import annotations
import json
from typing import Union, Any, Literal
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
from cherita.plotting.resampling import resample

MAX_SAMPLES = 100000
N_SAMPLES = 100000


def violin(
    adata_group: zarr.Group,
    mode: Literal["groupby", "multikey"],
    scale: Literal["width", "count"] = "width",
    **kwargs,
) -> Any:
    """Generate a Plotly violin plot JSON as a Python object from an
    Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object.
        mode (Literal["groupby", "multikey"]): Mode of the violin plot.
        scale (Literal["width", "count"], optional): Method to scale each violin's
            width.
            Can be set to "width" or "count". Defaults to "width".
        **kwargs: Additional keyword arguments.

    Returns:
        Any: A Plotly violin plot JSON as a Python object.
    """
    if mode == "groupby":
        fig, resampled = groupby_violin(adata_group, scale=scale, **kwargs)

    else:
        fig, resampled = multikey_violin(adata_group, scale=scale, **kwargs)

    fig.update_traces(scalemode=scale)
    fig_json = json.loads(fig.to_json())
    fig_json.update({"resampled": resampled})
    return fig_json


def groupby_violin(
    adata_group: zarr.Group,
    var_key: Union[str, dict],
    obs_col: str,
    obs_values: list[str] = None,
    scale: str = "width",
    var_names_col: str = None,
) -> tuple[go.Figure, bool]:
    """Generate a violin plot grouped by a categorical observation column.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object.
        var_key (Union[str, dict]): Marker key.
        obs_col (str): Name of the observation column to group by.
        obs_values (list[str], optional): List of observation values to include.
            Defaults to None.
        scale (str, optional): Method to scale each violin's width.
            Can be set to "width" or "count". Defaults to "width".
        var_names_col (str, optional): Name of the column in `adata_group.var`
            to use as violin names. Defaults to None.

    Returns:
        tuple[go.Figure, bool]: A tuple containing the Plotly violin plot figure
            and a boolean indicating if resampling was performed.
    """
    marker = Marker.from_any(adata_group, var_key, var_names_col)
    obs_colname = obs_col["name"]

    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    df = pd.DataFrame(
        marker.X,
        columns=[marker.name],
    )

    df[obs_colname], bins = to_categorical(obs, **obs_col)

    if obs_values is not None:
        df = df[df[obs_colname].isin(obs_values)]
        df[obs_colname] = df[obs_colname].cat.remove_unused_categories()

    violins = []
    resampled = False
    points = "outliers"
    for c in df[obs_colname].cat.categories:
        obs_data = df[marker.name][df[obs_colname] == c]
        if len(obs_data) >= MAX_SAMPLES:
            data_values = list(resample(obs_data, N_SAMPLES))
            resampled = True
            points = False
        else:
            data_values = obs_data

        violin = go.Violin(
            y=data_values,
            name=c,
            scalegroup="group" if scale == "count" else None,
            points=points,
        )
        violins.append(violin)

    fig = go.Figure(
        data=violins,
        layout=dict(
            yaxis=dict(title=marker.name),
            xaxis=dict(
                title=obs_colname + (f" ({bins} bins)" if bins else ""),
            ),
        ),
    )

    return fig, resampled


def multikey_violin(
    adata_group: zarr.Group,
    var_keys: list[Union[str, dict]] = [],
    obs_keys: list[str] = [],
    scale: str = "width",
    var_names_col: str = None,
) -> tuple[go.Figure, bool]:
    """Generate a multi-key violin plot.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object.
        var_keys (list[Union[str, dict]]): Key(s) to use for the violin plot.
        obs_keys (list[str]): List of observation keys.
        scale (str, optional): Method to scale each violin's width.
            Can be set to "width" or "count". Defaults to "width".
        var_names_col (str, optional): Name of the column in `adata_group.var`
            to use as violin names. Defaults to None.

    Returns:
        tuple[go.Figure, bool]: A tuple containing the Plotly violin plot figure
            and a boolean indicating if resampling was performed.
    """

    markers = [Marker.from_any(adata_group, m, var_names_col) for m in var_keys]

    if len(markers):
        df = pd.DataFrame(
            np.concatenate([[m.X] for m in markers]).T,
            columns=[m.name for m in markers],
        )
    else:
        df = pd.DataFrame()

    # Only numerical obs
    for o in obs_keys:
        if not isinstance(adata_group.obs[o], zarr.Array) and adata_group.obs[
            o
        ].dtype in [
            "int",
            "float",
        ]:
            raise BadRequest("Obs column '{}' is not numerical".format(o))
        df[o] = adata_group.obs[o]

    if not len(df.columns):
        raise BadRequest("No var or numerical obs provided")

    violins = []
    resampled = False
    points = "outliers"
    for col in df.columns:
        if len(df[col]) >= MAX_SAMPLES:
            data_values = list(resample(df[col], N_SAMPLES))
            resampled = True
            points = False
        else:
            data_values = df[col]

        violin = go.Violin(
            name=col,
            y=data_values,
            scalegroup="group" if scale == "count" else None,
            points=points,
        )
        violins.append(violin)

    fig = go.Figure(data=violins, layout=dict(yaxis=dict(title="Value")))

    return fig, resampled
