import re
from typing import Union
import zarr
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from cherita.utils.adata_utils import (
    get_group_index_name,
    parse_data,
    type_category,
    type_discrete,
    type_numeric,
    type_bool,
    continuous2categorical,
)
from cherita.plotting.resampling import resample

parse_dtype = {
    "category": type_category,
    "object": type_discrete,
    "int": type_numeric,
    "float": type_numeric,
    "complex": type_numeric,
    "bool": type_bool,
}

COL_NAME = "name"
INDEX_NAME = "matrix_index"


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


def get_obs_bin_data(
    adata_group: zarr.Group, obs_col: str, thresholds: list, nBins: int
):
    obs_s = parse_data(adata_group.obs[obs_col])
    cat = continuous2categorical(obs_s, thresholds, nBins)
    return type_category(pd.Series(cat))


def get_var_col_names(adata_group: zarr.Group):
    var_col_names = adata_group.var.attrs["column-order"]
    return var_col_names


# @TODO: optional obs_categories input
def get_var_histograms(
    adata_group: zarr.Group,
    var_index: int,
    obs_indices: Union[list, dict] = None,
):
    hist, bin_edges = np.histogram(adata_group.X[obs_indices or slice(None), var_index])
    bin_edges = [
        [float(bin_edges[i]), float(bin_edges[i + 1])]
        for i in range(len(bin_edges) - 1)
    ]
    log10_hist = np.log10(hist + 1)
    return {"hist": hist.tolist(), "bin_edges": bin_edges, "log10": log10_hist.tolist()}


def get_var_names(adata_group: zarr.Group, col: str = None):
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


def match_var_names(
    adata_group: zarr.Group, data_df: pd.DataFrame, col: str = None, right_key=COL_NAME
):
    var_df = pd.DataFrame.from_records(get_var_names(adata_group, col))
    matched_df = var_df.merge(
        data_df, how="right", left_on=COL_NAME, right_on=right_key
    )
    matched_df["index"].fillna(matched_df["gene_id"], inplace=True)
    matched_df[COL_NAME].fillna(matched_df["gene_name"], inplace=True)
    matched_df[INDEX_NAME].fillna(-1, inplace=True)
    return matched_df


def get_obsm_keys(adata_group: zarr.Group):
    obsm_keys = list(adata_group.obsm.keys())
    return obsm_keys


def get_kde_values(data):
    x = np.linspace(data.min(), data.max(), 1000)
    kde = gaussian_kde(data)
    kde_values = kde.evaluate(x)
    return x, kde_values


def get_obs_distribution(adata_group: zarr.Group, obs_colname: str):
    MAX_SAMPLES = 100000
    N_SAMPLES = 100000

    obs = parse_data(adata_group.obs[obs_colname])

    resampled = False
    if len(obs) >= MAX_SAMPLES:
        data_values = np.array(resample(obs, N_SAMPLES))
        resampled = True
    else:
        data_values = obs

    log_data_values = np.log(np.square(data_values))

    kde_values = get_kde_values(data_values)
    log_kde_values = get_kde_values(log_data_values)

    return {
        "kde_values": [kde_values[0].tolist(), kde_values[1].tolist()],
        "log_kde_values": [log_kde_values[0].tolist(), log_kde_values[1].tolist()],
        "resampled": resampled,
    }
