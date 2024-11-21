from flask import request, jsonify, Response
from flask_restx import Resource, fields, Namespace
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
    pseudospatial_masks,
)

ns = Namespace("plotting", description="Plotting operations", path="/")

heatmap_model = ns.model(
    "HeatmapModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "varKeys": fields.List(
            fields.String, required=True, description="List of selected markers"
        ),
        "obsCol": fields.Raw(required=True, description="Selected obs column"),
        "obsValues": fields.List(
            fields.String, required=False, description="Selected obs values"
        ),
        "varNamesCol": fields.String(description="Var names column"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/heatmap")
class Heatmap(Resource):
    @ns.doc(
        description=(
            "Generate a heatmap from the given AnnData object"
            "for the selected markers and obs column"
        ),
        responses={
            200: "Success",
            400: "Bad request",
            500: "Internal server error",
        },
    )
    @ns.expect(heatmap_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                heatmap(
                    adata_group=adata_group,
                    var_keys=json_data["varKeys"],
                    obs_col=json_data["obsCol"],
                    obs_values=json_data.get("obsValues"),
                    var_names_col=json_data.get("varNamesCol"),
                    obs_indices=json_data.get("obsIndices"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


dotplot_model = ns.model(
    "DotplotModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "varKeys": fields.List(
            fields.String, required=True, description="List of selected markers"
        ),
        "obsCol": fields.Raw(required=True, description="Selected obs column"),
        "obsValues": fields.List(
            fields.String, required=False, description="Selected obs values"
        ),
        "meanOnlyExpressed": fields.Boolean(description="Mean only expressed"),
        "expressionCutoff": fields.Float(description="Expression cutoff"),
        "standardScale": fields.String(description="Standard scale"),
        "varNamesCol": fields.String(description="Var names column"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/dotplot")
class Dotplot(Resource):
    @ns.doc(
        description=(
            "Generate a dotplot from the given AnnData object"
            "for the selected markers and obs column"
        ),
        responses={200: "Success", 400: "Bad request", 500: "Internal server error"},
    )
    @ns.expect(dotplot_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                dotplot(
                    adata_group=adata_group,
                    var_keys=json_data["varKeys"],
                    obs_col=json_data["obsCol"],
                    obs_values=json_data.get("obsValues"),
                    mean_only_expressed=json_data.get("meanOnlyExpressed", False),
                    expression_cutoff=json_data.get("expressionCutoff", 0.0),
                    standard_scale=json_data.get("standardScale"),
                    var_names_col=json_data.get("varNamesCol"),
                    obs_indices=json_data.get("obsIndices"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


matrixplot_model = ns.model(
    "MatrixplotModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "varKeys": fields.List(
            fields.String, required=True, description="List of selected markers"
        ),
        "obsCol": fields.Raw(required=True, description="Selected obs column"),
        "obsValues": fields.List(
            fields.String, required=False, description="Selected obs values"
        ),
        "standardScale": fields.String(description="Standard scale"),
        "varNamesCol": fields.String(description="Var names column"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/matrixplot")
class Matrixplot(Resource):
    @ns.doc(
        description=(
            "Generate a matrixplot from the given AnnData object"
            "for the selected markers and obs column"
        ),
        responses={200: "Success", 400: "Bad request", 500: "Internal server error"},
    )
    @ns.expect(matrixplot_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(
                matrixplot(
                    adata_group=adata_group,
                    var_keys=json_data["varKeys"],
                    obs_col=json_data["obsCol"],
                    obs_values=json_data.get("obsValues"),
                    standard_scale=json_data.get("standardScale"),
                    var_names_col=json_data.get("varNamesCol"),
                    obs_indices=json_data.get("obsIndices"),
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


violin_model = ns.model(
    "ViolinModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "mode": fields.String(required=True, description="Mode"),
        "scale": fields.String(description="Scale"),
        "varNamesCol": fields.String(description="Var names column"),
    },
)

multikey_violin_model = ns.inherit(
    "MultikeyViolinModel",
    violin_model,
    {
        "varKeys": fields.List(fields.String, description="List of var keys"),
        "obsKeys": fields.List(fields.String, description="List of obs keys"),
    },
)

groupby_violin_model = ns.inherit(
    "GroupbyViolinModel",
    violin_model,
    {
        "varKey": fields.String(description="Var key"),
        "obsCol": fields.String(description="Obs column"),
        "obsValues": fields.List(fields.String, description="List of obs values"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/violin")
class Violin(Resource):
    @ns.doc(
        description=(
            "Generate a violin plot from the given AnnData object"
            "for the selected keys and obs column"
        ),
        responses={200: "Success", 400: "Bad request", 500: "Internal server error"},
    )
    @ns.expect(violin_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            mode = json_data["mode"]
            if mode == "multikey":
                ns.expect(multikey_violin_model)
                params = dict(
                    var_keys=json_data.get("varKeys", []),
                    obs_keys=json_data.get("obsKeys", []),
                )
            elif mode == "groupby":
                ns.expect(groupby_violin_model)
                params = dict(
                    var_key=json_data["varKey"],
                    obs_col=json_data["obsCol"],
                    obs_values=json_data.get("obsValues"),
                    obs_indices=json_data.get("obsIndices"),
                )
            return jsonify(
                violin(
                    adata_group=adata_group,
                    mode=mode,
                    scale=json_data.get("scale", "width"),
                    var_names_col=json_data.get("varNamesCol"),
                    **params,
                )
            )
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


pseudospatial_gene_model = ns.model(
    "PseudospatialGeneModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "varKey": fields.String(
            required=True,
            description=(
                "Var key to plot, e.g. gene name or gene symbol."
                " The key must exist in the AnnData.var index"
                "or in the `varNamesCol` column if provided"
            ),
        ),
        "format": fields.String(
            description="Plot format. Can be 'png', 'html', or 'json'. Defaults to `png`"
        ),
        "maskSet": fields.String(
            description=(
                "Mask set to use. Must exist in AnnData.uns.masks."
                " Defaults to `spatial`"
            )
        ),
        "maskValues": fields.List(
            fields.String,
            description="Mask values to plot. If left unset, all masks will be plotted",
        ),
        "varNamesCol": fields.String(
            description="Column in var that contains `varKey` names"
        ),
        "obsCol": fields.Raw(
            description="Optional obs column to filter data by. To be used alongside `obsValues`"
        ),
        "obsValues": fields.List(
            fields.String, description="Obs values to filter data by"
        ),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
        "colormap": fields.String(description="Colormap name"),
        "fullHtml": fields.Boolean(
            description=(
                "Whether to return a full HTML document" " or just the `div` element"
            )
        ),
        "showColorbar": fields.Boolean(
            description="Whether to include a colorbar in the plot"
        ),
        "minValue": fields.Float(description="Min value for the color scale"),
        "maxValue": fields.Float(description="Max value for the color scale"),
        "width": fields.Integer(description="Plot width in pixels. Defaults to 500"),
        "height": fields.Integer(description="Plot height in pixels. Defaults to 500"),
    },
)


@ns.route("/pseudospatial/gene")
class PseudospatialGene(Resource):
    @ns.doc(
        description=(
            "Generate a pseudospatial plot for the given gene in the AnnData object"
        ),
        responses={
            200: "Success",
            400: "Bad request",
            500: "Internal server error",
            404: "Not found",
        },
    )
    @ns.expect(pseudospatial_gene_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            var_key = json_data["varKey"]
            plot_format = json_data.get("format", "png")
            optional_params = dict(
                mask_set="maskSet",
                mask_values="maskValues",
                var_names_col="varNamesCol",
                obs_col="obsCol",
                obs_values="obsValues",
                colormap="colormap",
                full_html="fullHtml",
                show_colorbar="showColorbar",
                min_value="minValue",
                max_value="maxValue",
                width="width",
                height="height",
                obs_indices="obsIndices",
            )
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_gene(
                adata_group,
                var_key=var_key,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            elif plot_format == "json":
                return jsonify(plot)
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


pseudospatial_categorical_model = ns.model(
    "PseudospatialCategoricalModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "obsCol": fields.Raw(required=True, description="Obs column to plot"),
        "obsValues": fields.List(
            fields.String,
            required=False,
            description="Obs values to plot. If left unset, all values will be plotted",
        ),
        "format": fields.String(description="Plot format"),
        "maskSet": fields.String(description="Mask set"),
        "mode": fields.String(
            description=(
                "Plotting mode. `across` for percentage across masks"
                "or `within` for percentage within masks"
            )
        ),
        "maskValues": fields.List(fields.String, description="Mask values"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
        "colormap": fields.String(description="Colormap"),
        "fullHtml": fields.Boolean(description="Full HTML"),
        "showColorbar": fields.Boolean(description="Show colorbar"),
        "minValue": fields.Float(description="Min value"),
        "maxValue": fields.Float(description="Max value"),
        "width": fields.Integer(description="Width"),
        "height": fields.Integer(description="Height"),
    },
)


@ns.route("/pseudospatial/categorical")
class PseudospatialCategorical(Resource):
    @ns.doc(
        description=(
            "Generate a pseudospatial plot for the given categorical obs column"
            "in the AnnData object"
        ),
        responses={
            200: "Success",
            400: "Bad request",
            500: "Internal server error",
            404: "Not found",
        },
    )
    @ns.expect(pseudospatial_categorical_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_col = json_data["obsCol"]
            plot_format = json_data.get("format", "png")

            optional_params = dict(
                obs_values="obsValues",
                mask_set="maskSet",
                mode="mode",
                mask_values="maskValues",
                colormap="colormap",
                full_html="fullHtml",
                show_colorbar="showColorbar",
                min_value="minValue",
                max_value="maxValue",
                width="width",
                height="height",
                obs_indices="obsIndices",
            )
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_categorical(
                adata_group,
                obs_col,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            elif plot_format == "json":
                return jsonify(plot)
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


pseudospatial_continuous_model = ns.model(
    "PseudospatialContinuousModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "obsCol": fields.Raw(required=True, description="Obs column"),
        "obsValues": fields.List(
            fields.String,
            required=False,
            description="Obs values to plot. If left unset, all values will be plotted",
        ),
        "format": fields.String(description="Plot format"),
        "maskSet": fields.String(description="Mask set"),
        "maskValues": fields.List(fields.String, description="Mask values"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
        "colormap": fields.String(description="Colormap"),
        "fullHtml": fields.Boolean(description="Full HTML"),
        "showColorbar": fields.Boolean(description="Show colorbar"),
        "minValue": fields.Float(description="Min value"),
        "maxValue": fields.Float(description="Max value"),
        "width": fields.Integer(description="Width"),
        "height": fields.Integer(description="Height"),
    },
)


@ns.route("/pseudospatial/continuous")
class PseudospatialContinuous(Resource):
    @ns.doc(
        description=(
            "Generate a pseudospatial plot for the given continuous obs column"
            "in the AnnData object"
        ),
        responses={
            200: "Success",
            400: "Bad request",
            500: "Internal server error",
            404: "Not found",
        },
    )
    @ns.expect(pseudospatial_continuous_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_col = json_data["obsCol"]
            plot_format = json_data.get("format", "png")

            optional_params = dict(
                obs_values="obsValues",
                mask_set="maskSet",
                mask_values="maskValues",
                colormap="colormap",
                full_html="fullHtml",
                show_colorbar="showColorbar",
                min_value="minValue",
                max_value="maxValue",
                width="width",
                height="height",
                obs_indices="obsIndices",
            )
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_continuous(
                adata_group,
                obs_col,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            elif plot_format == "json":
                return jsonify(plot)
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


pseudospatial_masks_model = ns.model(
    "PseudospatialMasksModel",
    {
        "url": fields.String(required=True, description="URL to the AnnData-Zarr file"),
        "format": fields.String(description="Plot format"),
        "maskSet": fields.String(description="Mask set"),
        "maskValues": fields.List(fields.String, description="Mask values"),
        "fullHtml": fields.Boolean(description="Full HTML"),
        "width": fields.Integer(description="Width"),
        "height": fields.Integer(description="Height"),
    },
)


@ns.route("/pseudospatial/masks")
class PseudospatialMasks(Resource):
    @ns.doc(
        description=(
            "Generate a pseudospatial plot of only the masks in the AnnData object"
        ),
        responses={
            200: "Success",
            400: "Bad request",
            500: "Internal server error",
            404: "Not found",
        },
    )
    @ns.expect(pseudospatial_masks_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            plot_format = json_data.get("format", "png")

            optional_params = {
                "mask_set": "maskSet",
                "full_html": "fullHtml",
                "width": "width",
                "height": "height",
            }
            optional_params_dict = {
                p: json_data.get(n)
                for p, n in optional_params.items()
                if n in json_data
            }

            plot = pseudospatial_masks(
                adata_group,
                plot_format=plot_format,
                **optional_params_dict,
            )

            if plot_format == "html":
                return Response(
                    plot,
                    mimetype="text/html",
                )
            elif plot_format == "json":
                return jsonify(plot)
            return plot
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")
