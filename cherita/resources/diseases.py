from flask import request, jsonify
from flask_restx import Resource, fields, Namespace
from cherita.resources.errors import BadRequest

from cherita.extensions import cache
from cherita.utils.adata_utils import open_anndata_zarr
from cherita.utils.caching import make_cache_key
from cherita.dataset.search import (
    search_diseases,
    search_disease_genes,
    get_disease_genes,
    get_disease_gene,
)

ns = Namespace("disease data", description="Disease related data", path="/")

disease_search_model = ns.model(
    "DiseaseSearchModel",
    {
        "diseaseDatasets": fields.List(
            fields.String,
            description="List of disease datasets to search in",
            required=True,
        ),
        "text": fields.String(description="Text to search for in disease names"),
    },
)


@ns.route("/diseases")
class SearchDiseaseNames(Resource):
    @ns.doc(
        description="Search for diseases in the given datasets",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(disease_search_model)
    def post(self):
        json_data = request.get_json()
        try:
            text = json_data.get("text", "")
            disease_datasets = json_data["diseaseDatasets"]
            return jsonify(search_diseases(disease_datasets, text))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


disease_genes_model = ns.model(
    "DiseaseGenesModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "diseaseName": fields.String(description="Name of the disease", required=True),
        "col": fields.String(description="Column to search for disease name"),
        "diseaseDatasets": fields.List(
            fields.String,
            description="List of disease datasets to search in",
            required=True,
        ),
    },
)


@ns.route("/disease/genes")
class GetDiseaseGenes(Resource):
    @ns.doc(
        description="Get genes associated with the given disease",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(disease_genes_model)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            disease_id = json_data["diseaseId"]
            col = json_data.get("col", None)
            disease_datasets = json_data.get("diseaseDatasets", [])
            return jsonify(
                get_disease_genes(adata_group, disease_id, col, disease_datasets)
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


disease_gene_model = ns.model(
    "DiseaseGeneModel",
    {
        "geneName": fields.String(description="Name of the gene", required=True),
        "diseaseDatasets": fields.List(
            fields.String,
            description="List of disease datasets to search in",
            required=True,
        ),
    },
)


@ns.route("/disease/gene")
class GetDiseaseGene(Resource):
    @ns.doc(
        description="Get diseases associated with the given gene",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(disease_gene_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            gene_name = json_data["geneName"]
            disease_datasets = json_data.get("diseaseDatasets", [])
            return jsonify(get_disease_gene(gene_name, disease_datasets))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


disease_gene_search_model = ns.model(
    "DiseaseGeneSearchModel",
    {
        "url": fields.String(description="URL to the zarr file", required=True),
        "text": fields.String(description="Text to search for in gene names"),
        "col": fields.String(description="Column to search for gene name"),
        "diseaseDatasets": fields.List(
            fields.String,
            description="List of disease datasets to search in",
            required=True,
        ),
    },
)


@ns.route("/diseases/genes")
class SearchDiseaseGenes(Resource):
    @ns.doc(
        description="Search for genes associated with diseases in the given datasets",
        responses={200: "Success", 400: "Invalid input", 500: "Internal server error"},
    )
    @ns.expect(disease_gene_search_model)
    @cache.cached(make_cache_key=make_cache_key, timeout=3600 * 24 * 7)
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            text = json_data.get("text", "")
            col = json_data.get("col", None)
            disease_datasets = json_data.get("diseaseDatasets", "")
            return jsonify(
                search_disease_genes(adata_group, col, text, disease_datasets)
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))
