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
