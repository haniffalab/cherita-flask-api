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
                )
            )
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


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
                    mean_only_expressed=json_data["meanOnlyExpressed"]
                    if "meanOnlyExpressed" in json_data
                    else False,
                    expression_cutoff=json_data["expressionCutoff"]
                    if "expressionCutoff" in json_data
                    else 0.0,
                    standard_scale=json_data["standardScale"]
                    if "standardScale" in json_data
                    else None,
                )
            )
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


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
                    standard_scale=json_data["standardScale"]
                    if "standardScale" in json_data
                    else None,
                )
            )
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class Violin(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                violin(
                    adata_group=adata_group,
                    keys=json_data["keys"],
                    obs_col=json_data["selectedObs"]
                    if "selectedObs" in json_data
                    else None,
                    scale=json_data["scale"] if "scale" in json_data else None,
                )
            )
        except BadRequest as e:
            raise e
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
