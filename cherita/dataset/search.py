import zarr
import pandas as pd
from cherita.dataset.metadata import get_var_names, match_var_names, COL_NAME
from cherita.utils.strapi_utils import get_from_strapi


def search_var_names(adata_group: zarr.Group, col: str = None, text: str = ""):
    var_df = pd.DataFrame.from_records(get_var_names(adata_group, col))
    var_df = var_df[var_df[COL_NAME].str.lower().str.startswith(text.lower())]

    return var_df.to_dict("records", index=True)


def search_diseases(
    text: str = "",
    disease_datasets: list = [],
):
    ENDPOINT = "disease-genes"

    params = {
        "filters[disease_datasets][name][$contains]": disease_datasets,
        f"filters[disease_name][${'startsWith' if len(text) < 2 else 'contains'}]": text,
        "fields": ["disease_id", "disease_name", "uid"],
        "pagination[start]": 0,
        "pagination[limit]": 500,
    }

    res = get_from_strapi(ENDPOINT, params)
    if "error" in res and res["error"]:
        return {"error": "Error fetching diseases"}

    res_data = [({"id": item["id"]} | item["attributes"]) for item in res["data"]]

    return res_data


# @TODO: change disease_name for disease_id
def get_disease_genes(
    adata_group: zarr.Group,
    disease_name: str,
    col: str = None,
    disease_datasets: list = [],
):
    ENDPOINT = "disease-genes"

    params = {
        "filters[disease_datasets][name][$contains]": disease_datasets,
        "filters[disease_name][$eq]": disease_name,
        "fields": ["disease_id", "disease_name", "gene_id", "gene_name", "uid"],
        "pagination[start]": 0,
        "pagination[limit]": 500,
    }

    res = get_from_strapi(ENDPOINT, params)
    if "error" in res and res["error"]:
        return {"error": "Error fetching diseases"}

    res_data = {item["id"]: item["attributes"] for item in res["data"]}
    res_df = pd.DataFrame.from_dict(res_data, orient="index")

    matched_df = match_var_names(adata_group, res_df, col, right_key="gene_name")

    return matched_df.to_dict("records", index=True)


def search_disease_genes(
    adata_group: zarr.Group,
    col: str = None,
    text: str = "",
    disease_datasets: list = [],
):
    ENDPOINT = "disease-genes"

    params = {
        "filters[disease_datasets][name][$contains]": disease_datasets,
        "filters[gene_name][$startsWith]": text,
        "fields": ["gene_id", "gene_name", "disease_id", "disease_name", "uid"],
        "pagination[start]": 0,
        "pagination[limit]": 500,
    }

    res = get_from_strapi(ENDPOINT, params)
    if "error" in res and res["error"]:
        return {"error": "Error fetching diseases"}

    res_data = {item["id"]: item["attributes"] for item in res["data"]}
    res_df = pd.DataFrame.from_dict(res_data, orient="index")

    matched_df = match_var_names(adata_group, res_df, col, right_key="gene_name")

    return matched_df.to_dict("records", index=True)
