from __future__ import annotations
from typing import Union
import zarr
from cherita.utils.models import Marker


# @TODO: optional obs_categories input
def get_var_x_mean(
    adata_group: zarr.Group,
    var_keys: list[Union[int, str, dict]],
    obs_indices: list[int] = None,
    var_names_col: str = None,
):
    markers = [
        Marker.from_any(adata_group, v, var_names_col=var_names_col) for v in var_keys
    ]
    return {m.name: float(m.get_X_at(obs_indices).mean()) for m in markers}
