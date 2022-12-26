"""
Settings for the development environment
"""

# Load secret key, database, etc from .env file
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import, wrong-import-position

# Development settings
DEBUG = True

EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
