import os


class Config(object):
    SECRET_KEY = os.environ.get("SECRET_KEY") or "best-password-here"
    API_VERSION = os.environ.get("API_VERSION") or "0.1.0"
    API_PREFIX = os.environ.get("API_PREFIX") or "/api"
    STRAPI_API = os.environ.get("STRAPI_API") or "http://localhost:1337/api/"
    DOCS_ROUTE = os.environ.get("DOCS_ROUTE") or "/docs/"
