from flask import Blueprint
from flask_restful import Api

from cherita.resources.about import About
from cherita.resources.errors import errors
from cherita.resources.plotting import (
    Heatmap,
    Dotplot,
    Matrixplot,
    Violin,
    PseudospatialGene,
    PseudospatialCategorical,
    PseudospatialContinuous,
)
from cherita.resources.dataset import (
    ObsColsNames,
    ObsCols,
    ObsBinData,
    ObsDistribution,
    ObsmKeys,
    VarColsNames,
    VarNames,
    VarHistograms,
)
from cherita.resources.diseases import (
    SearchDiseaseNames,
    SearchDiseaseGenes,
    GetDiseaseGenes,
    GetDiseaseGene,
)

bp = Blueprint("api_v1", __name__)
api = Api(bp, errors=errors)

api.add_resource(About, "/about")
api.add_resource(ObsColsNames, "/obs/cols/names")
api.add_resource(ObsCols, "/obs/cols")
api.add_resource(ObsBinData, "/obs/bins")
api.add_resource(ObsDistribution, "/obs/distribution")
api.add_resource(ObsmKeys, "/obsm/keys")
api.add_resource(VarColsNames, "/var/cols/names")
api.add_resource(VarNames, "/var/names")
api.add_resource(Heatmap, "/heatmap")
api.add_resource(Dotplot, "/dotplot")
api.add_resource(Matrixplot, "/matrixplot")
api.add_resource(Violin, "/violin")
api.add_resource(PseudospatialGene, "/pseudospatial/gene")
api.add_resource(PseudospatialCategorical, "/pseudospatial/categorical")
api.add_resource(PseudospatialContinuous, "/pseudospatial/continuous")
api.add_resource(SearchDiseaseNames, "/diseases")
api.add_resource(SearchDiseaseGenes, "/diseases/genes")
api.add_resource(GetDiseaseGenes, "/disease/genes")
api.add_resource(GetDiseaseGene, "/disease/gene")
api.add_resource(VarHistograms, "/var/histograms")
