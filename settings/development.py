# Load secret key, database, etc from .env file
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())

# Base settings
from .base import *

# Development settings
DEBUG = True

ALLOWED_HOSTS = []

# Email
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
EMAIL_HOST = os.environ['EMAIL_HOST']
EMAIL_PORT = os.environ['EMAIL_PORT']
EMAIL_HOST_USER = os.environ['EMAIL_HOST_USER']
EMAIL_HOST_PASSWORD = os.environ['EMAIL_HOST_PASSWORD']
EMAIL_USE_TLS = True
DEFAULT_FROM_EMAIL = os.environ['EMAIL_FROM']
