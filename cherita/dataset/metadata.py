import zarr
import pandas as pd
from cherita.utils.adata_utils import get_group_index_name, parse_data


def get_obs_col_names(adata_group: zarr.Group):
    obs_col_names = adata_group.obs.attrs["column-order"]
    return obs_col_names


def get_var_col_names(adata_group: zarr.Group):
    var_col_names = adata_group.var.attrs["column-order"]
    return var_col_names


def get_var_names(adata_group: zarr.Group, col: str = None):
    COL_NAME = "name"
    INDEX_NAME = "matrix_index"
    col = col or get_group_index_name(adata_group)
    var_df = pd.DataFrame(parse_data(adata_group.var[col]), columns=[COL_NAME])
    var_df.reset_index(names=[INDEX_NAME], inplace=True)
    var_df.sort_values(by=[COL_NAME], inplace=True)
    return var_df.to_dict("records", index=True)
