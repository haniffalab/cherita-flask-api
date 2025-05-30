from flask import request, jsonify
from flask_restx import Resource, fields, Namespace
from cherita.resources.errors import BadRequest

from cherita.extensions import cache
from cherita.utils.adata_utils import open_anndata_zarr
from cherita.utils.caching import make_cache_key
from cherita.dataset.metadata import (
    get_obs_col_histograms,
    get_obs_col_names,
    get_obs_col_metadata,
    get_pseudospatial_masks,
    get_var_col_names,
    get_obsm_keys,
    get_var_histograms,
    get_obs_bin_data,
    get_obs_distribution,
)
from cherita.dataset.matrix import get_var_x_mean
from cherita.dataset.search import search_var_names

ns = Namespace("dataset", description="Dataset related data", path="/")

obs_cols_name_model = ns.model(
    "ObsColsNamesModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
    },
)


@ns.route("/obs/cols/names")
class ObsColsNames(Resource):
    @ns.doc(
        description="Get the names of the observation columns in the dataset",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obs_cols_name_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obs_col_names(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


obs_cols_model = ns.model(
    "ObsColsModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "obsParams": fields.Raw(
            description=(
                "Parameters for the observation columns. Where obsParams is a "
                "dictionary with the obs column name as the key and the parameters as "
                "the value.  Possible values are `bins` dictionary with "
                "`nBins` or `thresholds` to bin continuous data. "
            )
        ),
        "retbins": fields.Boolean(
            description=(
                "Whether to return the bins for continuous data. "
                "Default is True, which will return the bins."
            ),
            default=True,
        ),
    },
)


@ns.route("/obs/cols")
class ObsCols(Resource):
    @ns.doc(
        description="Get the metadata of the observation columns in the dataset",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obs_cols_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_params = json_data.get("obsParams", {})
            retbins = json_data.get("retbins", True)
            return jsonify(
                get_obs_col_metadata(
                    adata_group, obs_params=obs_params, retbins=retbins
                )
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


obs_bin_data_model = ns.model(
    "ObsBinDataModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "obsCol": fields.String(
            description="Name of the observation column", required=True
        ),
        "thresholds": fields.List(fields.Float, description="Thresholds for binning"),
        "nBins": fields.Integer(description="Number of bins"),
    },
)


@ns.route("/obs/bins")
class ObsBinData(Resource):
    @ns.doc(
        description="Get the binned data for the observation column",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obs_bin_data_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_col = json_data["obsCol"]
            thresholds = json_data.get("thresholds", None)
            nBins = json_data.get("nBins", 5)
            if not thresholds and not nBins:
                raise BadRequest("Missing required parameter: thresholds or nBins")
            return jsonify(get_obs_bin_data(adata_group, obs_col, thresholds, nBins))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


obs_distribution_model = ns.model(
    "ObsDistributionModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "obsColname": fields.String(
            description="Name of the observation column", required=True
        ),
    },
)


@ns.route("/obs/distribution")
class ObsDistribution(Resource):
    @ns.doc(
        description="Get the distribution of the observation column",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obs_distribution_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            obs_colname = json_data["obsColname"]
            return jsonify(get_obs_distribution(adata_group, obs_colname))
        except KeyError as e:
            raise BadRequest(f"Missing required parameter: {e}")


obsm_keys_model = ns.model(
    "ObsmKeysModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
    },
)


@ns.route("/obsm/keys")
class ObsmKeys(Resource):
    @ns.doc(
        description="Get the keys of the observation metadata in the dataset",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obsm_keys_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_obsm_keys(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


var_cols_name_model = ns.model(
    "VarColsNamesModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
    },
)


@ns.route("/var/cols/names")
class VarColsNames(Resource):
    @ns.doc(
        description="Get the names of the variable columns in the dataset",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(var_cols_name_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            return jsonify(get_var_col_names(adata_group))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


var_names_model = ns.model(
    "VarNamesModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "col": fields.String(description="Name of the variable column"),
        "text": fields.String(description="Text to search for in variable names"),
    },
)


@ns.route("/var/names")
class VarNames(Resource):
    @ns.doc(
        description="Search for variable names in the dataset",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(var_names_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            col = json_data.get("col", None)
            text = json_data.get("text", "")
            return jsonify(search_var_names(adata_group, col, text))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


var_histograms_model = ns.model(
    "VarHistogramsModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "varKey": fields.Integer(description="Index of the variable"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/var/histograms")
class VarHistograms(Resource):
    @ns.doc(
        description="Get the histograms of the variable",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(var_histograms_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            var_key = json_data["varKey"]
            obs_indices = json_data.get("obsIndices", None)
            return jsonify(get_var_histograms(adata_group, var_key, obs_indices))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


obs_histograms_model = ns.model(
    "ObsHistogramsModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "varKey": fields.Integer(description="Index of the variable"),
        "obsCol": fields.Raw(required=True, description="Selected obs column"),
        "obsIndices": fields.List(
            fields.Integer, description="List of observation indices"
        ),
    },
)


@ns.route("/obs/histograms")
class ObsHistograms(Resource):
    @ns.doc(
        description="Get the histograms of the observation column",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(obs_histograms_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            var_key = json_data["varKey"]
            obs_col = json_data["obsCol"]
            obs_indices = json_data.get("obsIndices", None)
            return jsonify(
                get_obs_col_histograms(adata_group, var_key, obs_col, obs_indices)
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


# @TODO: add optional filters by obs or indices
matrix_mean_model = ns.model(
    "VarXMeanModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "varKeys": fields.List(fields.String, description="List of variable keys"),
        "varNamesCol": fields.String(description="Name of the variable names column"),
    },
)


@ns.route("/matrix/mean")
class MatrixMean(Resource):
    @ns.doc(
        description="Get the mean of the variables",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(matrix_mean_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            var_keys = json_data["varKeys"]
            var_names_col = json_data.get("varNamesCol", None)
            obs_indices = json_data.get("obsIndices", None)
            return jsonify(
                get_var_x_mean(
                    adata_group,
                    var_keys,
                    obs_indices=obs_indices,
                    var_names_col=var_names_col,
                )
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


@ns.route("/masks")
class Masks(Resource):
    @ns.doc(
        description="Get the pseudospatial masks data",
        responses={200: "Success", 404: "Not found", 500: "Internal server error"},
    )
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        adata_group = open_anndata_zarr(json_data["url"])
        return jsonify(get_pseudospatial_masks(adata_group))
