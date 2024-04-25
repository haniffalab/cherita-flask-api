from __future__ import annotations
from flask import current_app
import logging
import requests
from urllib.parse import urljoin


def get_from_strapi(endpoint: str, params: dict = {}):

    API_URL = current_app.config["STRAPI_API"]

    response = requests.get(urljoin(API_URL, endpoint), params=params)
    if not response:
        logging.error("Error fetching disease genes")
        return {"error": "Error fetching disease genes"}
    if response.status_code != 200:
        logging.warning(
            f"Unsuccesful status code {response.status_code} fetching disease genes"
        )
        return {"error": response}

    return response.json()
