# Base settings
from .base import *

# Test settings
DEBUG = False

ALLOWED_HOSTS = [os.environ['ALLOWED_HOST']]

# Email
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
EMAIL_HOST = 'example.com'
EMAIL_PORT = 25
EMAIL_HOST_USER = 'user'
EMAIL_HOST_PASSWORD = 'password'
EMAIL_USE_TLS = True
DEFAULT_FROM_EMAIL = 'user@example.com'