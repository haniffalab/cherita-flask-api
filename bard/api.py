from flask import Blueprint
from flask_restful import Api

from bard.resources.about import About
from bard.resources.errors import errors
from bard.resources.plotting import Heatmap
from bard.resources.dataset import ObsCols, VarCols, VarNames

bp = Blueprint("api_v1", __name__)
api = Api(bp, errors=errors)

api.add_resource(About, "/about")
api.add_resource(Heatmap, "/heatmap")
api.add_resource(ObsCols, "/obs/cols")
api.add_resource(VarCols, "/var/cols")
api.add_resource(VarNames, "/var/names")
