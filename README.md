[![dependabot](https://img.shields.io/badge/dependabot-enabled-025E8C?logo=dependabot&logoColor=white)](https://github.com/haniffalab/cherita-flask-api/security)
[![python](https://shields.io/badge/python-3.11-blue)](https://github.com/haniffalab/cherita-flask-api?tab=readme-ov-file#development)

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

### Caching

The API is setup to use Redis for caching (supported by [Flask-Caching](https://flask-caching.readthedocs.io/en/latest/)).
The app will attempt to connect to a Redis instance if the environment variable `REDIS_HOST` is set to the instance's IP address. `REDIS_PORT` and `CACHE_KEY_PREFIX` can also be set.
If `REDIS_HOST` is set and the app cannot connect to the Redis instance it will throw an "Redis connection error" on each attempt to access the cache.

When updating the API's responses you will have to flush the cache to avoid getting outdated data.
You need to connect to the redis instance, we recommend using [redis-cli](https://redis.io/docs/latest/develop/tools/cli/).
You can install it with

```sh
sudo apt-get install redis-tools
```

You can send the `FLUSHALL` command directly like

```sh
redis-cli -h instance-ip-address -p port FLUSHALL
```

Alternatively, you can connect to the instance and then execute the command

```sh
redis-cli -h instance-ip-address -p port
FLUSHALL
```

#### GCP Memorystore

Note that when using a [Memorystore Redis instance](https://cloud.google.com/memorystore/docs/redis/memorystore-for-redis-overview) you will need to connect from a VM that is within the instance's authorized network. Refer to the [official documentation](https://cloud.google.com/memorystore/docs/redis/connect-redis-instance#connecting-compute-engine-redis-cli) for more information.