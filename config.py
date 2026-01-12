import os

version = "1.2.11"


class Config(object):
    SECRET_KEY = os.environ.get("SECRET_KEY", "best-password-here")
    API_VERSION = os.environ.get("API_VERSION", version)
    API_PREFIX = os.environ.get("API_PREFIX", "/api")
    STRAPI_API = os.environ.get("STRAPI_API", "http://localhost:1337/api/")
    DOCS_ROUTE = os.environ.get("DOCS_ROUTE", "/docs/")
    DOCS_JSON = os.environ.get("DOCS_JSON", "docs/swagger.json")
    REDIS_HOST = os.environ.get("REDIS_HOST", None)
    REDIS_PORT = os.environ.get("REDIS_PORT", 6379)
    CACHE_KEY_PREFIX = os.environ.get("CACHE_KEY_PREFIX", "cherita-flask-cache_")
