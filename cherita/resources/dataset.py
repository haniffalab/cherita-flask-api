from flask import request, jsonify
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.dataset.metadata import (
    get_obs_col_names,
    get_obs_col_metadata,
    get_var_col_names,
    get_obsm_keys,
    get_var_histograms,
    get_obs_bin_data,
    get_obs_distribution,
)
from cherita.dataset.search import search_var_names


class ObsColsNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obs_col_names(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class ObsCols(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obs_col_metadata(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class ObsBinData(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_col = json_data["obs_col"]
            thresholds = json_data.get("thresholds", None)
            nBins = json_data.get("nBins", None)
            if not thresholds and not nBins:
                raise BadRequest("Missing required parameter: thresholds or nBins")
            return jsonify(get_obs_bin_data(adata_group, obs_col, thresholds, nBins))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class ObsDistribution(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_colname = json_data["obs_colname"]
            return jsonify(get_obs_distribution(adata_group, obs_colname))
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class ObsmKeys(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obsm_keys(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class VarColsNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_var_col_names(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class VarNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            col = json_data.get("col", None)
            text = json_data.get("text", "")
            return jsonify(search_var_names(adata_group, col, text))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class VarHistograms(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            var_index = json_data["var_index"]
            obs_indices = json_data.get("obs_indices", None)
            return jsonify(get_var_histograms(adata_group, var_index, obs_indices))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
