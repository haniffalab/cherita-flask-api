from flask import Blueprint
from flask_restful import Api

from cherita.resources.about import About
from cherita.resources.errors import errors


bp = Blueprint("api_v1", __name__)
api = Api(bp, errors=errors)

api.add_resource(About, "/about")

