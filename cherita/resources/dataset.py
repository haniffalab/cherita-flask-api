from flask import request, jsonify
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.dataset.metadata import (
    get_obs_col_names,
    get_obs_col_metadata,
    get_var_col_names,
    get_obsm_keys,
)
from cherita.dataset.search import search_var_names


class ObsColsNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obs_col_names(adata_group))
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
        except Exception as e:
            raise e


class ObsCols(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obs_col_metadata(adata_group))
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
        except Exception as e:
            raise e


class ObsmKeys(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obsm_keys(adata_group))
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
        except Exception as e:
            raise e


class VarColsNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_var_col_names(adata_group))
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
        except Exception as e:
            raise e


class VarNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            col = json_data.get("col", None)
            text = json_data.get("text", "")
            return jsonify(search_var_names(adata_group, col, text))
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
        except Exception as e:
            raise e
