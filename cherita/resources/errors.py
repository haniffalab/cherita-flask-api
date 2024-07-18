from werkzeug.exceptions import HTTPException


class InternalServerError(HTTPException):
    pass


class ReadZarrError(HTTPException):
    pass


class InvalidKey(HTTPException):
    pass


class InvalidObs(HTTPException):
    pass


class InvalidVar(HTTPException):
    pass


class BadRequest(HTTPException):
    pass


class FetchError(HTTPException):
    pass


class NotInData(HTTPException):
    pass


errors = {
    "InternalServerError": {
        "name": "InternalServerError",
        "message": "Something went wrong",
        "status": 500,
    },
    "ReadZarrError": {"name": "ReadZarrError", "status": 400},
    "InvalidKey": {"name": "InvalidKey", "status": 400},
    "InvalidObs": {"name": "InvalidObs", "status": 400},
    "InvalidVar": {"name": "InvalidVar", "status": 400},
    "BadRequest": {"name": "BadRequest", "status": 400},
    "FetchError": {"name": "FetchError", "status": 500},
    "NotInData": {"name": "NotInData", "status": 404},
}
