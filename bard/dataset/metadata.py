import zarr
import anndata as ad
from bard.utils.adata_utils import get_group_index_name, parse_data


def get_obs_col_names(adata_group: zarr.Group):
    obs_col_names = adata_group.obs.attrs["column-order"]
    return obs_col_names


def get_var_col_names(adata_group: zarr.Group):
    var_col_names = adata_group.var.attrs["column-order"]
    return var_col_names


def get_var_names(adata_group: zarr.Group, col: str = None):
    col = col or get_group_index_name(adata_group)
    var_names = list(parse_data(adata_group.var[col]))
    var_names.sort()
    return var_names
