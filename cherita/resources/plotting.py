from flask import request, jsonify
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.plotting.heatmap import heatmap
from cherita.plotting.dotplot import dotplot
from cherita.plotting.matrixplot import matrixplot
from cherita.plotting.violin import violin


class Heatmap(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                heatmap(
                    adata_group=adata_group,
                    markers=json_data["selectedMultiVar"],
                    obs_col=json_data["selectedObs"],
                    obs_values=json_data.get("obsValues"),
                    var_names_col=json_data.get("varNamesCol"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class Dotplot(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                dotplot(
                    adata_group=adata_group,
                    markers=json_data["selectedMultiVar"],
                    obs_col=json_data["selectedObs"],
                    obs_values=json_data.get("obsValues"),
                    mean_only_expressed=json_data.get("meanOnlyExpressed", False),
                    expression_cutoff=json_data.get("expressionCutoff", 0.0),
                    standard_scale=json_data.get("standardScale"),
                    var_names_col=json_data.get("varNamesCol"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class Matrixplot(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                matrixplot(
                    adata_group=adata_group,
                    markers=json_data["selectedMultiVar"],
                    obs_col=json_data["selectedObs"],
                    obs_values=json_data.get("obsValues"),
                    standard_scale=json_data.get("standardScale"),
                    var_names_col=json_data.get("varNamesCol"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class Violin(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                violin(
                    adata_group=adata_group,
                    keys=json_data["keys"],
                    obs_col=json_data.get("selectedObs"),
                    obs_values=json_data.get("obsValues"),
                    scale=json_data.get("scale"),
                    var_names_col=json_data.get("varNamesCol"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")
