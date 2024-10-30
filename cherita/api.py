from flask import Blueprint
from flask_restx import Api

from config import Config
import cherita.resources.errors as errors
from cherita.resources.about import About
from cherita.resources.plotting import ns as plotting_namespace
from cherita.resources.dataset import (
    ns as dataset_namespace,
)
from cherita.resources.diseases import (
    ns as diseases_namespace,
)

bp = Blueprint("api_v1", __name__)
api = Api(bp, doc=Config.DOCS_ROUTE, version=Config.API_VERSION)

api.add_resource(About, "/about")

api.add_namespace(dataset_namespace)
api.add_namespace(plotting_namespace)
api.add_namespace(diseases_namespace)


@api.errorhandler(errors.InternalServerError)
def handle_internal_server_error(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.ReadZarrError)
def handle_read_zarr_error(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.InvalidKey)
def handle_invalid_key(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.InvalidObs)
def handle_invalid_obs(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.InvalidVar)
def handle_invalid_var(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.BadRequest)
def handle_bad_request(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.FetchError)
def handle_fetch_error(error):
    return {"message": error.message}, error.status_code


@api.errorhandler(errors.NotInData)
def handle_not_in_data(error):
    return {"message": error.message}, error.status_code
