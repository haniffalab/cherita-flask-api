from flask import request, jsonify, Response
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.plotting.heatmap import heatmap
from cherita.plotting.dotplot import dotplot
from cherita.plotting.matrixplot import matrixplot
from cherita.plotting.violin import violin
from cherita.plotting.pseudospatial import (
    pseudospatial_gene,
    pseudospatial_categorical,
    pseudospatial_continuous,
)


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
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            plot_format = json_data.get("format", "png")
            optional_params = {
                "marker_id": "varId",
                "marker_name": "varName",
                "mask": "mask",
                "mask_values": "maskValues",
                "var_names_col": "varNamesCol",
                "colormap": "colormap",
                "full_html": "fullHtml",
                "show_colorbar": "showColorbar",
                "min_value": "minValue",
                "max_value": "maxValue",
                "width": "width",
                "height": "height",
            }
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_gene(
                adata_group, plot_format=plot_format, **optional_params_dict
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class PseudospatialCategorical(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_colname = json_data["obsName"]
            obs_values = json_data["obsValues"]
            plot_format = json_data.get("format", "png")

            optional_params = {
                "mask": "mask",
                "mode": "mode",
                "mask_values": "maskValues",
                "colormap": "colormap",
                "full_html": "fullHtml",
                "show_colorbar": "showColorbar",
                "min_value": "minValue",
                "max_value": "maxValue",
                "width": "width",
                "height": "height",
            }
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_categorical(
                adata_group,
                obs_colname,
                obs_values,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


class PseudospatialContinuous(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_colname = json_data["obsName"]
            plot_format = json_data.get("format", "png")

            optional_params = {
                "mask": "mask",
                "mask_values": "maskValues",
                "colormap": "colormap",
                "full_html": "fullHtml",
                "show_colorbar": "showColorbar",
                "min_value": "minValue",
                "max_value": "maxValue",
                "width": "width",
                "height": "height",
            }
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_continuous(
                adata_group,
                obs_colname,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")
