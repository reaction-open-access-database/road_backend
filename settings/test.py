"""
Settings for the testing environment
"""

import os

os.environ = {
    "EMAIL_HOST": "example.com",
    "EMAIL_PORT": "25",
    "EMAIL_HOST_USER": "user",
    "EMAIL_HOST_PASSWORD": "password",
    "EMAIL_FROM": "user@example.com",
    "ADMIN_URL": "admin/",
    "LOGIN_URL": "login/",
    "ALLOWED_HOST": "localhost",
} | os.environ

# If this is a development environment, load secret key, database, etc from .env file
from dotenv import find_dotenv, load_dotenv  # pylint: disable=wrong-import-position

load_dotenv(find_dotenv())

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import, wrong-import-position

# Test settings
DEBUG = False

# Email
EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
