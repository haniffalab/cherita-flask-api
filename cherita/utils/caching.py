import hashlib
import json

from flask import current_app, jsonify, request
from google.appengine.api import memcache


def generate_cache_key(data):
    data_str = json.dumps(data, sort_keys=True)
    return hashlib.md5(data_str.encode("utf-8")).hexdigest()


def cached(expiration=3600):
    """Decorator to cache the response of a function."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            if not current_app.config["IS_PRODUCTION"]:
                return func(*args, **kwargs)

            json_data = request.get_json()
            cache_key = generate_cache_key(json_data)

            cached_result = memcache.get(cache_key)
            if cached_result:
                return jsonify(cached_result)

            jsonified_result = func(*args, **kwargs)

            memcache.set(cache_key, json.loads(jsonified_result), expiration)

            return jsonified_result

        return wrapper

    return decorator
