import hashlib
import json

from flask import current_app, jsonify, request
from google.appengine.api import memcache


def generate_cache_key():
    data = {
        "method": request.method,
        "path": request.path,
        "body": request.get_json(silent=True) or {},
    }
    data_str = json.dumps(data, sort_keys=True)
    return hashlib.md5(data_str.encode("utf-8")).hexdigest()


def cached(expiration=3600):
    """Decorator to cache the response of a function."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            if not current_app.config["IS_PRODUCTION"]:
                return func(*args, **kwargs)

            cache_key = generate_cache_key()

            cached_result = memcache.get(cache_key)
            if cached_result:
                return jsonify(json.loads(cached_result))

            response = func(*args, **kwargs)

            memcache.set(cache_key, json.dumps(response.get_json()), expiration)

            return response

        return wrapper

    return decorator
