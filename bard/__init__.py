import os

from flask import Flask, render_template
from werkzeug.middleware.proxy_fix import ProxyFix

from config import Config


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

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    @app.route("/")
    def index():
        return render_template("index.html")

    # apply the blueprints to the app
    from bard import api
    app.register_blueprint(
        api.bp,
        url_prefix=app.config["API_PREFIX"]
    )

    return app
