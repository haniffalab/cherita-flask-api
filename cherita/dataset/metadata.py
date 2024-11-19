import logging
import re
from typing import Union
import zarr
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from cherita.resources.errors import InvalidObs, NotInData
from cherita.utils.adata_utils import (
    get_group_index_name,
    parse_data,
    type_category,
    type_discrete,
    type_numeric,
    type_bool,
    to_categorical,
)
from cherita.plotting.resampling import resample
from cherita.utils.models import Marker

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
    cat, _ = to_categorical(
        obs_s, type="continuous", thresholds=thresholds, nBins=nBins, fillna=False
    )
    return {**type_category(pd.Series(cat)), "type": "continuous"}


def get_var_col_names(adata_group: zarr.Group):
    var_col_names = adata_group.var.attrs["column-order"]
    return var_col_names


def histogram(a: np.ndarray, hist_range: tuple[float, float] = None):
    hist, bin_edges = np.histogram(a, range=hist_range)
    bin_edges = [
        [float(bin_edges[i]), float(bin_edges[i + 1])]
        for i in range(len(bin_edges) - 1)
    ]
    log10_hist = np.log10(hist + 1)
    return {"hist": hist.tolist(), "bin_edges": bin_edges, "log10": log10_hist.tolist()}


def get_var_histograms(
    adata_group: zarr.Group,
    var_key: Union[int, str, dict],
    obs_indices: list[int] = None,
):
    marker = Marker.from_any(adata_group, var_key)
    return histogram(marker.get_X_at(obs_indices))


def get_obs_col_histograms(
    adata_group: zarr.Group,
    var_key: Union[int, str, dict],
    obs_col: dict,
    obs_indices: list[int] = None,
):
    obs_colname = obs_col["name"]
    if obs_colname not in adata_group.obs:
        raise NotInData(f"Column '{obs_colname}' not found in AnnData")
    try:
        obs = parse_data(adata_group.obs[obs_colname])
    except KeyError as e:
        raise InvalidObs(f"Invalid observation {e}")

    categorical_obs, _ = to_categorical(obs, **obs_col)
    categorical_obs = (
        categorical_obs[obs_indices] if obs_indices is not None else categorical_obs
    )

    marker = Marker.from_any(adata_group, var_key)
    X = marker.get_X_at(obs_indices)
    min_X, max_X = X.min(), X.max()

    if min_X == max_X:
        max_X += 1

    obs_data = {}
    for cat in categorical_obs.categories:
        obs_data[cat] = histogram(
            X[np.flatnonzero(categorical_obs == cat)], [min_X, max_X]
        )

    return obs_data


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
    data = data[~np.isnan(data)]
    data = data[~np.isinf(data)]

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

    log_data_values = np.log(np.square(data_values) + np.spacing(1))

    kde_values = get_kde_values(data_values)
    log_kde_values = get_kde_values(log_data_values)

    return {
        "kde_values": [kde_values[0].tolist(), kde_values[1].tolist()],
        "log_kde_values": [log_kde_values[0].tolist(), log_kde_values[1].tolist()],
        "resampled": resampled,
    }


def get_pseudospatial_masks(adata_group: zarr.Group):
    if "masks" not in adata_group.uns:
        raise NotInData("masks not found in adata_group.uns")

    mask_data = {}
    for mask in adata_group.uns["masks"]:
        m = adata_group.uns["masks"][mask]
        if "polygons" not in m or "obs" not in m:
            logging.err(f"polygons not found in adata_group.uns['masks'][{mask}]")
        else:
            obs = parse_data(adata_group.obs[m["obs"][()]])
            mask_data[mask] = obs.categories.tolist()

    if not mask_data:
        raise NotInData("No masks found in adata_group.uns['masks']")

    return mask_data
