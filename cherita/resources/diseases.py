from flask import request, jsonify
from flask_restful import Resource
from cherita.resources.errors import BadRequest

from cherita.utils.adata_utils import open_anndata_zarr
from cherita.dataset.search import (
    search_diseases,
    search_disease_genes,
    get_disease_genes,
    get_disease_gene,
)


class SearchDiseaseNames(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            text = json_data.get("text", "")
            disease_datasets = json_data.get("diseaseDatasets", [])
            return jsonify(search_diseases(text, disease_datasets))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class GetDiseaseGenes(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            disease_name = json_data["diseaseName"]
            col = json_data.get("col", None)
            disease_datasets = json_data.get("diseaseDatasets", [])
            return jsonify(
                get_disease_genes(adata_group, disease_name, col, disease_datasets)
            )
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class GetDiseaseGene(Resource):
    def post(self):
        json_data = request.get_json()
        try:
            adata_group = open_anndata_zarr(json_data["url"])
            gene_name = json_data["geneName"]
            disease_datasets = json_data.get("diseaseDatasets", [])
            return jsonify(get_disease_gene(adata_group, gene_name, disease_datasets))
        except KeyError as e:
            raise BadRequest("Missing required parameter: {}".format(e))


class SearchDiseaseGenes(Resource):
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
