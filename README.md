[![python-tests](https://github.com/haniffalab/webatlas-pipeline/actions/workflows/tests-python.yml/badge.svg)](https://github.com/haniffalab/webatlas-pipeline/actions/workflows/tests-python.yml)
[![codecov](https://codecov.io/gh/haniffalab/webatlas-pipeline/branch/main/graph/badge.svg?token=7HQVFH08WJ)](https://app.codecov.io/gh/haniffalab/webatlas-pipeline)
[![dependabot](https://img.shields.io/badge/dependabot-enabled-025E8C?logo=dependabot&logoColor=white)](https://github.com/haniffalab/cherita-flask-api/security)
[![python](https://shields.io/badge/python-3.10-blue)](https://github.com/haniffalab/cherita-flask-api?tab=readme-ov-file#development)

# Cherita Flask API

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://haniffalab.github.io/webatlas-pipeline)
[![demo](https://img.shields.io/badge/Demos-view-blue)](https://cellatlas.io/webatlas)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.7405818.svg)](https://doi.org/10.5281/zenodo.7405818)

Desc

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
