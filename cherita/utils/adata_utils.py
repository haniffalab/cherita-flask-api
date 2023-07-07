from __future__ import annotations
import os
import zarr
import numpy as np
import pandas as pd
from zarr.errors import GroupNotFoundError
from typing import Union
from urllib.parse import urlparse, ParseResult

from cherita.resources.errors import BadRequest


def get_group_index(group: zarr.Group):
    if "_index" in group.attrs:
        return group[group.attrs["_index"]]
    else:
        return group["_index"]


def get_group_index_name(group: zarr.Group):
    if "_index" in group.attrs:
        return group.attrs["_index"]
    else:
        return "_index"


def get_index_in_array(array: zarr.Array, item: str):
    return np.where(array[:] == item)[0][0]


def get_indices_in_array(array: zarr.Array, items: list[str]):
    sorter = np.argsort(array[:])
    return sorter[np.searchsorted(array[:], items, sorter=sorter)]


def parse_data(data: Union[zarr.Group, zarr.Array], store: zarr.Group = None):
    if type(data) == zarr.Group:
        return parse_group(data)
    elif type(data) == zarr.Array:
        return parse_array(data, store)


def parse_group(group: zarr.Group):
    if "_index" in group.attrs:
        df = pd.DataFrame(index=group[group.attrs["_index"]])
        for name in [
            name
            for name in group.array_keys()
            if not name.startswith("_") and name != get_group_index_name(group)
        ]:
            df[name] = group[name]
        for name in group.group_keys():
            df[name] = parse_group(group[name])
        return df
    elif "codes" in group and "categories" in group:
        series = pd.Categorical.from_codes(
            group["codes"][:], categories=group["categories"][:]
        )
        if np.array_equal(
            np.sort(series.categories.values), np.sort(["True", "False"])
        ):
            return series.map({"True": True, "False": False}).astype(bool)
        return series


def parse_array(array: zarr.Array, store: zarr.Group = None):
    if "categories" in array.attrs:
        series = pd.Categorical.from_codes(
            array[:],
            categories=store[
                os.path.join(os.path.dirname(array.path), array.attrs["categories"])
            ],
        )
        if np.array_equal(
            np.sort(series.categories.values), np.sort(["True", "False"])
        ):
            return series.map({"True": True, "False": False}).astype(bool)
        return series
    else:
        return array[:]


def get_s3_http_options(o: ParseResult):
    bucket, endpoint = o.hostname.split(".", 1)
    s3url = "s3://" + bucket + o.path
    storage_options = {"anon": True, "endpoint_url": o.scheme + "://" + endpoint}
    return s3url, storage_options


def open_anndata_zarr(url: str):
    o = urlparse(url)

    if o.scheme in ["gcs", "gs"]:
        storage_options = {"token": "anon"}
    elif o.netloc == "storage.googleapis.com":
        storage_options = None
    else:
        url, storage_options = get_s3_http_options(o)

    try:
        adata_group = zarr.open_consolidated(
            url, storage_options=storage_options, mode="r"
        )
    except (FileNotFoundError, KeyError):
        try:
            adata_group = zarr.open_group(
                url, storage_options=storage_options, mode="r"
            )
        except GroupNotFoundError:
            raise BadRequest("Cannot open Anndata Zarr at URL {}".format(url))

    return adata_group


def type_category(obs):
    categories = [str(i) for i in obs.cat.categories.values.flatten()]

    if len(categories) > 100:
        return {
            "type": "categorical",
            "is_truncated": True,
            "values": categories[:99],
        }

    return {
        "type": "categorical",
        "is_truncated": False,
        "values": categories,
    }


def type_bool(obs):
    return {"type": "categorical", "values": {1: "True", 0: "False"}}


def type_numeric(obs):
    accuracy = 4
    return {
        "type": "continuous",
        "min": encode_dtype(round(ndarray_min(obs), accuracy)),
        "max": encode_dtype(round(ndarray_max(obs), accuracy)),
        "mean": encode_dtype(round(ndarray_mean(obs), accuracy)),
        "median": encode_dtype(round(ndarray_median(obs), accuracy)),
    }


def type_discrete(obs):
    return {"type": "discrete"}


def ndarray_max(a):
    if np.isnan(a).all():
        return 0
    else:
        return np.nanmax(a)


def ndarray_min(a):
    if np.isnan(a).all():
        return 0
    else:
        return np.nanmin(a)


def ndarray_mean(a):
    if np.isnan(a).all():
        return 0
    else:
        return np.nanmean(a)


def ndarray_median(a):
    if np.isnan(a).all():
        return 0
    else:
        return np.nanmedian(a)


def encode_dtype(a):
    if hasattr(a, "dtype"):
        if isinstance(a, np.integer):
            return int(a)
        if isinstance(a, np.floating):
            return float(a)
        if isinstance(a, np.bool_):
            return bool(a)
        if isinstance(a, np.ndarray) or isinstance(a, pd.Categorical):
            return a.tolist()
        return a
    else:
        return a
