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


def get_obs_col_metadata(
    adata_group: zarr.Group,
    cols: [str] = None,
    obs_params: dict = {},
    retbins: bool = True,
):
    obs_df = parse_data(adata_group.obs)
    obs_metadata = []

    cols = cols or obs_df.columns

    for col in cols:
        if col not in obs_df.columns:
            continue
        t = re.sub(r"[^a-zA-Z]", "", obs_df[col].dtype.name)
        metadata = parse_dtype[t](
            obs_df[col], obs_params=obs_params.get(col, {}), retbins=retbins
        )
        if t == "category":
            colors = get_obs_colors(adata_group, col)
            if colors and metadata["n_values"] == len(colors):
                metadata["colors"] = colors

        obs_metadata.append(
            {
                "name": col,
                **metadata,
            }
        )
    return obs_metadata


def get_obs_bin_data(
    adata_group: zarr.Group, obs_col: str, thresholds: list, nBins: int = 5
):
    obs_s = parse_data(adata_group.obs[obs_col])
    cat, _ = to_categorical(
        obs_s, type="continuous", thresholds=thresholds, nBins=nBins, fillna=False
    )
    return {**type_category(pd.Series(cat)), "type": "continuous"}


def is_hex_triplet(color: str) -> bool:
    if not isinstance(color, str):
        return False
    return bool(re.fullmatch(r"^#[0-9a-fA-F]{6}$", color))


def get_obs_colors(adata_group: zarr.Group, obs_col: str) -> list[str] | None:
    # Following cxg schema https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/6.0.0/schema.md#column_colors
    obs_colors = f"{obs_col}_colors"
    if obs_colors in adata_group.uns.array_keys():
        try:
            colors = parse_data(adata_group.uns[obs_colors]).tolist()
            if not all(is_hex_triplet(c) for c in colors):
                logging.warn(f"Colors in {obs_colors} are not hex strings")
                return None
            return colors

        except Exception as e:
            logging.error(f"Failed to parse colors for {obs_col}: {e}")
            return None
    return None


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
    if len(X.shape) > 1:
        X = X.mean(0)
    if not X.shape[-1]:
        min_X, max_X = 0, 1
    else:
        min_X, max_X = X.min(), X.max()

    if min_X == max_X:
        max_X += 1

    obs_data = {}
    for cat in categorical_obs.categories:
        obs_data[cat] = histogram(
            X[np.flatnonzero(categorical_obs == cat)],
            [min_X, max_X],
        )

    return obs_data


def get_var_names(adata_group: zarr.Group, col: str = None, names: [str] = None):
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
    if names is not None:
        var_df = var_df[var_df[COL_NAME].isin(names)]
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


def get_obsm_keys(adata_group: zarr.Group, filter2d: bool = True):
    if filter2d:
        # Filter to only obsm keys that can be plotted in 2D
        # i.e. is a zarr array, has float, int or uint dtype, and has at least 2 dimensions
        return [
            key
            for key in adata_group.obsm.keys()
            if isinstance(adata_group.obsm[key], zarr.Array)
            and adata_group.obsm[key].dtype.kind in "fiu"
            and adata_group.obsm[key].ndim == 2
            and adata_group.obsm[key].shape[1] > 1
        ]
    else:
        return list(adata_group.obsm.keys())


def get_kde_values(data):
    data = data[~np.isnan(data)]
    data = data[~np.isinf(data)]

    if not len(data):
        return np.array([]), np.array([])

    x = np.linspace(data.min(), data.max(), 250)
    kde = gaussian_kde(data)
    kde_values = kde.evaluate(x)
    return x, kde_values


def get_obs_distribution(adata_group: zarr.Group, obs_colname: str):
    MAX_SAMPLES = 25000
    N_SAMPLES = 25000

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



def get_obs_values(adata_group: zarr.Group, col: str, values: [str] = None):
    """
    Get unique values from an obs column in the AnnData object.

    Parameters
    ----------
    adata_group : zarr.Group
        AnnData Zarr store
    col : str
        Name of the obs column
    values : list of str, optional
        Restrict to this subset of values

    Returns
    -------
    dict
        Records of unique values with an index
    """
    import pandas as pd

    if col not in adata_group.obs:  # or check with adata_group["obs"].keys()
        raise KeyError(f"Column '{col}' not found in .obs")

    # Parse the column values
    obs_vals = pd.Series(parse_data(adata_group.obs[col])).astype(str)

    # Deduplicate & sort
    # obs_df = pd.DataFrame({ "value": sorted(obs_vals.unique()) })
    obs_df = pd.DataFrame(
        {COL_NAME: sorted(obs_vals.unique())}
    )

    # if values is not None:
    #     obs_df = obs_df[obs_df["value"].isin(values)]
    if values is not None:
        obs_df = obs_df[obs_df[COL_NAME].isin(values)]

    # Reset index like get_var_names
    obs_df.reset_index(names=["index"], inplace=True)
    obs_df.reset_index(names=[INDEX_NAME], inplace=True)

    return obs_df.to_dict("records", index=True)
