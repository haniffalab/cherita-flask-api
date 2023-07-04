from __future__ import annotations
import os
import zarr
import numpy as np
import pandas as pd
from typing import Union
from urllib.parse import urlparse, ParseResult


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
    except:
        adata_group = zarr.open_group(url, storage_options=storage_options, mode="r")
    return adata_group
