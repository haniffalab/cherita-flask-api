from flask import request, jsonify
from flask_restful import Resource

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.dataset.metadata import get_obs_col_names, get_var_col_names, get_var_names


class ObsCols(Resource):
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        return jsonify(get_obs_col_names(adata_group))


class VarCols(Resource):
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        return jsonify(get_var_col_names(adata_group))


class VarNames(Resource):
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        col = json_data["col"] if "col" in json_data else None
        return jsonify(get_var_names(adata_group, col))
