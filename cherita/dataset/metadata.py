import re
import zarr
import pandas as pd
from cherita.utils.adata_utils import (
    get_group_index_name,
    parse_data,
    type_category,
    type_discrete,
    type_numeric,
    type_bool,
)

parse_dtype = {
    "category": type_category,
    "object": type_discrete,
    "int": type_numeric,
    "float": type_numeric,
    "complex": type_numeric,
    "bool": type_bool,
}


def get_obs_col_names(adata_group: zarr.Group):
    obs_col_names = adata_group.obs.attrs["column-order"]
    return obs_col_names


def get_obs_col_metadata(adata_group: zarr.Group):
    obs_df = parse_data(adata_group.obs)
    obs_metadata = []
    for col in obs_df.columns:
        t = re.sub(r"[^a-zA-Z]", "", obs_df[col].dtype.name)
        obs_metadata.append({"name": col, **parse_dtype[t](obs_df[col])})
    return obs_metadata


def get_var_col_names(adata_group: zarr.Group):
    var_col_names = adata_group.var.attrs["column-order"]
    return var_col_names


def get_var_names(adata_group: zarr.Group, col: str = None):
    COL_NAME = "name"
    INDEX_NAME = "matrix_index"
    idx_col = get_group_index_name(adata_group.var)
    col = col or idx_col
    var_df = pd.DataFrame(
        parse_data(adata_group.var[col]),
        columns=[COL_NAME],
        index=parse_data(adata_group.var[idx_col]),
    )
    var_df.reset_index(names=["index"], inplace=True)
    var_df.reset_index(names=[INDEX_NAME], inplace=True)
    var_df.sort_values(by=[COL_NAME], inplace=True)
    return var_df.to_dict("records", index=True)


def get_obsm_keys(adata_group: zarr.Group):
    obsm_keys = list(adata_group.obsm.keys())
    return obsm_keys
