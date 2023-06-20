from flask import request, jsonify
from flask_restful import Resource, reqparse

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.plotting.heatmap import heatmap


class Heatmap(Resource):
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        return jsonify(
            heatmap(
                adata_group, json_data["selectedMultiVar"], json_data["selectedObs"]
            )
        )
