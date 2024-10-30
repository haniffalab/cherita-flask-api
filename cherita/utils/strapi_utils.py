from __future__ import annotations
from flask import current_app
import requests
from urllib.parse import urljoin

from cherita.resources.errors import FetchError


def get_from_strapi(endpoint: str, params: dict = {}):
    API_URL = current_app.config["STRAPI_API"]

    response = requests.get(urljoin(API_URL, endpoint), params=params)
    if not response:
        raise FetchError("Error fetching from strapi")
    if response.status_code != 200:
        raise FetchError(
            f"Unsuccesful status code {response.status_code} fetching from strapi"
        )

    return response.json()
