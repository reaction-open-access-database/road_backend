"""
Settings for the testing environment
"""

import os

# If this is a development environment, load secret key, database, etc from .env file
from dotenv import find_dotenv, load_dotenv

load_dotenv(find_dotenv())

DEFAULT_ENVIRON = {
    "EMAIL_HOST": "example.com",
    "EMAIL_PORT": "25",
    "EMAIL_HOST_USER": "user",
    "EMAIL_HOST_PASSWORD": "password",
    "EMAIL_FROM": "user@example.com",
    "ADMIN_URL": "admin/",
    "LOGIN_URL": "/login/",
    "ALLOWED_HOST": "localhost",
}
for key, value in DEFAULT_ENVIRON.items():
    if key not in os.environ:
        os.environ[key] = value

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import, wrong-import-position

# Test settings
DEBUG = False

ALLOW_REMOTE_DATABASE_FLUSH = True

LOGGING |= {
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "simple",
        }
    },
    "loggers": {
        "road": {
            "handlers": ["console"],
            "level": "DEBUG",
            "propagate": True,
        }
    },
}

# Email
if os.environ.get("SEND_EMAILS", "False") != "True":
    EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
