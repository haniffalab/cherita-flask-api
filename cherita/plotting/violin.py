from __future__ import annotations
import json
from typing import Union, Any
import zarr
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from cherita.resources.errors import BadRequest, InvalidKey, InvalidObs, InvalidVar

from cherita.utils.adata_utils import (
    to_categorical,
    get_group_index,
    get_indices_in_array,
    parse_data,
)

MAX_SAMPLES = 100000
N_SAMPLES = 100000


def violin(
    adata_group: zarr.Group,
    keys: Union[list[int], list[str]],
    obs_col: dict = None,
    scale: str = "width",
    var_names_col: str = None,
) -> Any:
    """Method to generate a Plotly violin plot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        keys (Union[list[int], list[str]]): Keys of .var_names or numerical obs columns.
        obs_col (dict, optional): Obs colum to group by. Defaults to None.
        standard_scale (str, optional): Method to scale each violin's width.
            Can be set to "width" or "count".
            Defaults to "width".
        var_names_col (str, optional): Column in var to pull markers' names from.
            Defaults to None.

    Returns:
        Any: A Plotly violin plot JSON as a Python object
    """
    if obs_col:
        if not isinstance(obs_col, dict):
            raise BadRequest("'selectedObs' must be an object")
        if not isinstance(keys, str):
            raise BadRequest(
                "'keys' parameter should be a single string item"
                "when grouping by an observation"
            )

        if not all(isinstance(x, int) for x in keys) and not all(
            isinstance(x, str) for x in keys
        ):
            raise InvalidVar(
                "List of features should be all of the same type str or int"
            )

        if isinstance(keys[0], str):
            try:
                marker_idx = get_indices_in_array(
                    get_group_index(adata_group.var), keys
                )
            except InvalidKey:
                raise InvalidVar(f"Invalid feature name {keys}")
        else:
            marker_idx = keys

        if var_names_col:
            keys = adata_group.var[var_names_col][marker_idx]
        else:
            keys = get_group_index(adata_group.var)[marker_idx]

        obs_colname = obs_col["name"]
        try:
            obs = parse_data(adata_group.obs[obs_colname])
        except KeyError as e:
            raise InvalidObs(f"Invalid observation {e}")

        df = pd.DataFrame(adata_group.X[:, marker_idx], columns=[keys])

        df[obs_colname], bins = to_categorical(obs, **obs_col)

        violins = []
        resampled = False
        points = "outliers"
        for c in df[obs_colname].cat.categories:
            obs_data = df[keys][df[obs_colname] == c]
            if len(obs_data) >= MAX_SAMPLES:
                data_values = kde_resample(obs_data, N_SAMPLES)
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
                yaxis=dict(title=keys),
                xaxis=dict(
                    title=obs_colname + (f" ({bins} bins)" if bins else ""),
                ),
            ),
        )

    else:
        if isinstance(keys, str):
            keys = [keys]

        if not all(isinstance(x, str) for x in keys):
            raise BadRequest("'keys' should be of type str for multikey violin plot")

        var_keys = list(parse_data(adata_group.var).index.intersection(keys))
        obs_keys = list(set(adata_group.obs.attrs["column-order"]).intersection(keys))

        try:
            marker_idx = get_indices_in_array(
                get_group_index(adata_group.var), var_keys
            )
        except InvalidKey:
            raise InvalidVar(f"Invalid features {var_keys}")

        if var_names_col:
            var_keys = adata_group.var[var_names_col][marker_idx]

        df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=var_keys)

        # Only numerical obs
        for k in obs_keys:
            if not isinstance(adata_group.obs[k], zarr.Array) and adata_group.obs[
                k
            ].dtype in [
                "int",
                "float",
            ]:
                raise BadRequest("Obs column '{}' is not numerical".format(k))
            df[k] = adata_group.obs[k]

        violins = []
        resampled = False
        points = "outliers"
        for col in df.columns:
            if len(df[col]) >= MAX_SAMPLES:
                data_values = kde_resample(df[col], N_SAMPLES)
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

    fig.update_traces(scalemode=scale)
    fig_json = json.loads(fig.to_json())
    fig_json.update({"resampled": resampled})
    return fig_json


def kde_resample(df: pd.DataFrame, nsamples: int):
    np.random.seed(nsamples)
    kde_values = [df.min(), df.max()]
    unq, ids = np.unique(df, return_inverse=True)
    all_ids = np.random.choice(ids, size=10**6, replace=True)
    ar = np.bincount(all_ids) / 10**6
    kde_values.extend(np.random.choice(a=unq, size=nsamples, p=ar))
    return list(kde_values)
