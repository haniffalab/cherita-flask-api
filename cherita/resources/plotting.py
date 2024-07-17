from flask import request, jsonify
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.plotting.heatmap import heatmap
from cherita.plotting.dotplot import dotplot
from cherita.plotting.matrixplot import matrixplot
from cherita.plotting.violin import violin
from cherita.plotting.pseudospatial import pseudospatial_gene


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
                    var_names_col=json_data.get("varNamesCol", None),
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
                    mean_only_expressed=json_data.get("meanOnlyExpressed", False),
                    expression_cutoff=json_data.get("expressionCutoff", 0.0),
                    standard_scale=json_data.get("standardScale", None),
                    var_names_col=json_data.get("varNamesCol", None),
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
                    standard_scale=json_data.get("standardScale", None),
                    var_names_col=json_data.get("varNamesCol", None),
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
                    obs_col=json_data.get("selectedObs", None),
                    scale=json_data.get("scale", None),
                    var_names_col=json_data.get("varNamesCol", None),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class PseudospatialGene(Resource):
    def get(self):
        params = request.args
        try:
            adata_group = open_anndata_zarr(params["url"])
            optional_params = {
                "marker_id": "varId",
                "marker_name": "varName",
                "mask": "mask",
                "var_names_col": "varNamesCol",
                "colormap": "colormap",
                "plot_format": "format",
                "full_html": "fullHtml",
                "show_colorbar": "showColorbar",
                "min_value": "minValue",
                "max_value": "maxValue",
                "width": "width",
                "height": "height",
            }
            optional_params_dict = {
                p: params.get(n) for p, n in optional_params.items() if n in params
            }

            return pseudospatial_gene(adata_group, **optional_params_dict)
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")
