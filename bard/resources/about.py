from flask import current_app
from flask_restful import Resource


class About(Resource):
    def get(self):
        return {
            "name": current_app.name, 
            "version": current_app.config["API_VERSION"]
        }
