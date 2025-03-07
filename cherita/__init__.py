import os
import click
import logging
import json
from flask import Flask, render_template
from flask.cli import with_appcontext
from werkzeug.middleware.proxy_fix import ProxyFix

from flask_cors import CORS

from config import Config
from cherita.api import api, bp
from cherita.extensions import cache


def create_app(test_config=None):
    # create and configure an instance of the Flask application
    app = Flask(
        __name__,
        instance_relative_config=True,
        static_url_path="",
        static_folder="static",
        template_folder="templates",
    )

    # rewrite headers for proxy deployments
    app.wsgi_app = ProxyFix(
        app.wsgi_app, x_for=1, x_proto=0, x_host=0, x_port=0, x_prefix=1
    )

    # load config
    app.config.from_object(Config)
    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.update(test_config)

    if app.config["REDIS_HOST"] is not None:
        logging.info("Using RedisCache")
        cache.init_app(
            app,
            config={
                "CACHE_TYPE": "cherita.utils.caching.SafeRedisCache",
                "CACHE_REDIS_HOST": app.config["REDIS_HOST"],
                "CACHE_REDIS_PORT": app.config["REDIS_PORT"],
                "CACHE_KEY_PREFIX": app.config["CACHE_KEY_PREFIX"],
            },
        )
    else:
        logging.warning("No REDIS_HOST provided, using NullCache")
        cache.init_app(app, config={"CACHE_TYPE": "NullCache"})

    CORS(app)

    @app.route("/")
    def index():
        return render_template("index.html")

    # apply the blueprints to the app
    app.register_blueprint(bp, url_prefix=app.config["API_PREFIX"])

    @click.command("docs")
    @with_appcontext
    def generate_docs():
        with app.test_request_context():
            os.makedirs(os.path.dirname(app.config["DOCS_JSON"]), exist_ok=True)
            with open(app.config["DOCS_JSON"], "w") as f:
                json.dump(api.__schema__, f, indent=4)
            click.echo("Generated docs")

    app.cli.add_command(generate_docs)

    return app
