from flask import request, jsonify
from flask_restful import Resource, reqparse

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.plotting.heatmap import heatmap
from cherita.plotting.dotplot import dotplot


class Heatmap(Resource):
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        return jsonify(
            heatmap(
                adata_group=adata_group,
                markers=json_data["selectedMultiVar"],
                obs_col=json_data["selectedObs"],
            )
        )


class Dotplot(Resource):
    def post(self):
        json_data = request.get_json()
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
            )
        )
