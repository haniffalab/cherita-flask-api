from werkzeug.exceptions import HTTPException


class InternalServerError(HTTPException):
    pass


class BadRequest(HTTPException):
    pass


errors = {
    "InternalServerError": {"message": "Something went wrong", "status": 500},
    "BadRequest": {"status": 400},
}
