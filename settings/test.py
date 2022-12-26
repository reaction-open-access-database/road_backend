"""
Settings for the testing environment
"""

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import

# Test settings
DEBUG = False

ALLOWED_HOSTS = [os.environ["ALLOWED_HOST"]]

ADMIN_URL = "admin/"

# Email
EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
EMAIL_HOST = "example.com"
EMAIL_PORT = 25
EMAIL_HOST_USER = "user"
EMAIL_HOST_PASSWORD = "password"
EMAIL_USE_TLS = True
DEFAULT_FROM_EMAIL = "user@example.com"
