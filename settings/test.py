"""
Settings for the testing environment
"""

import os

os.environ["EMAIL_HOST"] = "example.com"
os.environ["EMAIL_PORT"] = "25"
os.environ["EMAIL_HOST_USER"] = "user"
os.environ["EMAIL_HOST_PASSWORD"] = "password"
os.environ["EMAIL_FROM"] = "user@example.com"

os.environ["ADMIN_URL"] = "admin/"
os.environ["LOGIN_URL"] = "login/"

os.environ["ALLOWED_HOST"] = "localhost"

# If this is a development environment, load secret key, database, etc from .env file
from dotenv import find_dotenv, load_dotenv  # pylint: disable=wrong-import-position

load_dotenv(find_dotenv())

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import, wrong-import-position

# Test settings
DEBUG = False

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
EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
