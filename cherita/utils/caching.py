import hashlib
import json
import redis
import logging

from flask import request
from flask_caching.backends.rediscache import RedisCache


def make_cache_key(*args, **kwargs) -> str:
    data = {
        "method": request.method,
        "path": request.path,
        "body": request.get_json(silent=True) or {},
    }
    data_str = json.dumps(data, sort_keys=True)
    return hashlib.md5(data_str.encode("utf-8")).hexdigest()


class SafeRedisCache(RedisCache):
    def _log_connection_error(self, e):
        logging.error("Redis connection error", e)

    def get(self, *args, **kwargs):
        try:
            return super().get(*args, **kwargs)
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return None

    def set(self, *args, **kwargs):
        try:
            return super().set(*args, **kwargs)
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return False

    def delete(self, *args, **kwargs):
        try:
            return super().delete(*args, **kwargs)
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return False

    def clear(self):
        try:
            return super().clear()
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return False

    def inc(self):
        try:
            return super().inc()
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return False

    def dec(self):
        try:
            return super().dec()
        except redis.ConnectionError as e:
            self._log_connection_error(e)
            return False
