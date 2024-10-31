class InternalServerError(Exception):
    def __init__(self, message="Internal server error", status_code=500):
        self.message = message
        self.status_code = status_code


class ReadZarrError(Exception):
    def __init__(self, message="Error reading zarr", status_code=400):
        self.message = message
        self.status_code = status_code


class InvalidKey(Exception):
    def __init__(self, message="Invalid key", status_code=400):
        self.message = message
        self.status_code = status_code


class InvalidObs(Exception):
    def __init__(self, message="Invalid obs", status_code=400):
        self.message = message
        self.status_code = status_code


class InvalidVar(Exception):
    def __init__(self, message="Invalid var", status_code=400):
        self.message = message
        self.status_code = status_code


class BadRequest(Exception):
    def __init__(self, message="Bad request", status_code=400):
        self.message = message
        self.status_code = status_code


class FetchError(Exception):
    def __init__(self, message="Error fetching data", status_code=500):
        self.message = message
        self.status_code = status_code


class NotInData(Exception):
    def __init__(self, message="Data not found", status_code=404):
        self.message = message
        self.status_code = status_code
