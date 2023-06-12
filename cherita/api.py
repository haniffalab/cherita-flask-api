from flask import Blueprint
from flask_restful import Api

from cherita.resources.about import About
from cherita.resources.errors import errors
from cherita.resources.plotting import Heatmap
from cherita.resources.dataset import ObsCols, VarCols, VarNames

bp = Blueprint("api_v1", __name__)
api = Api(bp, errors=errors)

api.add_resource(About, "/about")
api.add_resource(Heatmap, "/heatmap")
api.add_resource(ObsCols, "/obs/cols")
api.add_resource(VarCols, "/var/cols")
api.add_resource(VarNames, "/var/names")
