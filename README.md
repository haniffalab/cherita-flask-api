## Development

Requires Python3.10

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
