import hashlib
import json

from flask import request


def make_cache_key(*args, **kwargs) -> str:
    data = {
        "method": request.method,
        "path": request.path,
        "body": request.get_json(silent=True) or {},
    }
    data_str = json.dumps(data, sort_keys=True)
    return hashlib.md5(data_str.encode("utf-8")).hexdigest()
