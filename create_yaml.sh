#!/bin/sh

FILE="app.yaml"
/bin/cat <<EOM >$FILE
runtime: python310
service: $APP_SERVICE
service_account: $SERVICE_ACCOUNT
instance_class: $INSTANCE_CLASS
entrypoint: gunicorn -b :$PORT -w 2 'cherita:create_app()'
env_variables:
  FLASK_APP: $FLASK_APP
automatic_scaling:
  min_instances: $MIN_INSTANCES
  max_instances: $MAX_INSTANCES
  target_cpu_utilization: $TARGET_CPU_UTILIZATION
  max_concurrent_requests: $MAX_CONCURRENT_REQUESTS
EOM
