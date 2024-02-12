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
    if isinstance(data, zarr.Group):
        return parse_group(data)
    elif isinstance(data, zarr.Array):
        return parse_array(data, store)


def parse_group(group: zarr.Group):
    encoding_type = group.attrs.get("encoding-type", "")
    if encoding_type == "dataframe":
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
    elif encoding_type == "categorical":
        if "codes" in group and "categories" in group:
            series = pd.Categorical.from_codes(
                group["codes"][:], categories=group["categories"][:]
            )
            if np.array_equal(
                np.sort(series.categories.values), np.sort(["True", "False"])
            ):
                return series.map({"True": True, "False": False}).astype(bool)
            return series
        else:
            raise SystemError(
                f"Categorical group {group} does not contain 'codes' and 'categories'"
            )
    elif encoding_type == "dict":
        d = {}
        for k in group.keys():
            d[k] = parse_data(group[k])
        return d
    else:
        raise SystemError(f"Unrecognized encoding-type {encoding_type}")


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

    return {
        "type": "categorical",
        "is_truncated": True,
        "values": categories[:99] if len(categories) > 100 else categories,
        "n_values": len(categories),
        "state": {},
    }


def type_bool(obs):
    return {"type": "bool", "values": {1: "True", 0: "False"}}


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
    obs = obs.astype("category")
    categories = [str(i) for i in obs.cat.categories.values.flatten()]

    return {
        "type": "discrete",
        "is_truncated": True,
        "values": categories[:99] if len(categories) > 100 else categories,
        "n_values": len(categories),
        "state": {},
    }


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


def to_categorical(data: Union[pd.Series, np.Array], type: str, **kwargs):
    func_dict = {
        "continuous": continuous2categorical,
        "discrete": discrete2categorical,
        "bool": bool2categorical,
    }

    return func_dict.get(type, lambda data, **kwargs: (data, None))(data, **kwargs)


def continuous2categorical(
    s: pd.Series,
    thresholds: list[Union[int, float]] = None,
    nBins: int = 20,
    start: int = 1,
    **kwargs,
):
    return (
        pd.Categorical(pd.cut(s, thresholds, include_lowest=True, labels=False) + start)
        if thresholds
        else pd.Categorical(pd.cut(s, nBins, include_lowest=True, labels=False))
    ), nBins


def discrete2categorical(array: np.Array, nBins: int = 20, start: int = 1, **kwargs):
    s = pd.Series(array).astype("category")
    if nBins >= len(s.cat.categories):
        return s, None
    else:
        return (
            pd.Categorical(
                pd.cut(s.index, nBins, include_lowest=True, labels=False) + start
            ),
            nBins,
        )


def bool2categorical(array: np.Array, **kwargs):
    return pd.Series(array).astype("category"), None
