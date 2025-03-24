from __future__ import annotations
from typing import Union
import zarr
import numpy as np
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
    X_dict = {m.name: m.get_X_at(obs_indices) for m in markers}
    return {k: float(np.mean(v)) if v.shape[-1] else None for k, v in X_dict.items()}
