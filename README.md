[![dependabot](https://img.shields.io/badge/dependabot-enabled-025E8C?logo=dependabot&logoColor=white)](https://github.com/haniffalab/cherita-flask-api/security)
[![python](https://shields.io/badge/python-3.10-blue)](https://github.com/haniffalab/cherita-flask-api?tab=readme-ov-file#development)

# Cherita Flask API

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://cherita-flask-api-dot-haniffa-lab.nw.r.appspot.com/api/docs/)
[![demo](https://img.shields.io/badge/Demo-view-blue)](https://cherita-flask-api-dot-haniffa-lab.nw.r.appspot.com/)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.14772563.svg)](https://doi.org/10.5281/zenodo.14772563)

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
