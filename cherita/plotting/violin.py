from __future__ import annotations
import json
from typing import Union, Any
import zarr
import pandas as pd
import plotly.graph_objects as go

from cherita.utils.adata_utils import (
    get_group_index,
    get_index_in_array,
    get_indices_in_array,
    parse_data,
)


def violinplot(
    adata_group: zarr.Group,
    keys: Union[str, list[str]],
    obs_col: str = None,
    scale: str = "width",
) -> Any:
    """Method to generate a Plotly violin plot JSON as a Python object
    from an Anndata-Zarr object.

    Args:
        adata_group (zarr.Group): Root zarr Group of an Anndata-Zarr object
        keys (list[str], optional): Keys of .var_names or numerical obs columns.
        obs_col (str, optional): Obs colum to group by. Defaults to None.
        standard_scale (str, optional): Method to scale each violin's width.
            Can be set to "width" or "count".
            Defaults to "width".

    Returns:
        Any: A Plotly violin plot JSON as a Python object
    """
    if isinstance(keys, str) and obs_col:
        marker_idx = get_index_in_array(get_group_index(adata_group.var), keys)
        df = pd.DataFrame(adata_group.X[:, marker_idx], columns=[keys])

        obs = parse_data(adata_group.obs[obs_col])
        df[obs_col] = obs

        violins = []
        for c in df[obs_col].cat.categories:
            violin = go.Violin(y=df[keys][df[obs_col] == c], name=c)
            violins.append(violin)

        fig = go.Figure(
            data=violins, layout=dict(yaxis=dict(title=keys), xaxis=dict(title=obs_col))
        )

    elif isinstance(keys, str):
        keys = [keys]

        var_keys = list(parse_data(adata_group.var).index.intersection(keys))
        obs_keys = list(set(adata_group.obs.attrs["column-order"]).intersection(keys))

        marker_idx = get_indices_in_array(get_group_index(adata_group.var), var_keys)
        df = pd.DataFrame(adata_group.X.oindex[:, marker_idx], columns=var_keys)

        # Only numerical obs
        # @TODO: return error
        for k in obs_keys:
            if isinstance(adata_group.obs[k], zarr.Array) and adata_group.obs[
                k
            ].dtype in [
                "int",
                "float",
            ]:
                df[k] = adata_group.obs[k]

        violins = []
        for col in df.columns:
            violin = go.Violin(name=col, y=df[col])
            violins.append(violin)

        fig = go.Figure(violins)

    return json.loads(fig.to_json())
