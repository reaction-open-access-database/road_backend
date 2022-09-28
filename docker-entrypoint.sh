#!/bin/bash

python manage.py collectstatic --noinput
# TODO: Make the migrations before committing
python manage.py makemigrations
python manage.py migrate
exec "$@"