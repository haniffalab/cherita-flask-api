from flask import current_app
from flask_restx import Resource


class About(Resource):
    def get(self):
        return {"name": current_app.name, "version": current_app.config["API_VERSION"]}
