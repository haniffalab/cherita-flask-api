[![dependabot](https://img.shields.io/badge/dependabot-enabled-025E8C?logo=dependabot&logoColor=white)](https://github.com/haniffalab/cherita-flask-api/security)
[![python](https://shields.io/badge/python-3.10-blue)](https://github.com/haniffalab/cherita-flask-api?tab=readme-ov-file#development)

# Cherita Flask API

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://cherita-flask-api-925164757806.europe-west2.run.app/api/docs/)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.14772563.svg)](https://doi.org/10.5281/zenodo.14772563)

A web service designed to provide a robust and scalable interface for managing and analyzing biological data. This API leverages Flask, a lightweight WSGI web application framework in Python, to deliver high-performance endpoints for data retrieval, processing, and visualization.

## Development

Create a development environment with the `requirements.txt` file

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Set the `FLASK_APP` environment variable to the `cherita` directory

```sh
export FLASK_APP=cherita
```

Run the API

```sh
python -m flask run
```

You can enable debugging mode with

```sh
python -m flask run --debug
```

You can set a specific port for the API with the `-p` flag

```sh
python -m flask run -p 8001
```
