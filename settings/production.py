# Base settings
from .base import *

# Production settings
DEBUG = False

ALLOWED_HOSTS = ['reaction-open-access-database.herokuapp.com']

# SSL
CSRF_COOKIE_SECURE = True
SESSION_COOKIE_SECURE = True
SECURE_SSL_REDIRECT = True

# HSTS
SECURE_HSTS_SECONDS = 60
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = True
