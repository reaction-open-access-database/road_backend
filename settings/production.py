"""
Settings for the production environment
"""

# Base settings
from .base import *  # pylint: disable=wildcard-import, unused-wildcard-import

# Production settings
DEBUG = False

SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")

LOGGING |= {
    "handlers": {
        "file": {
            "level": "DEBUG",
            "class": "logging.handlers.RotatingFileHandler",
            "filename": LOG_DIR / "django.log",
            "maxBytes": 16 * 1024 * 1024,  # 16 MB
            "formatter": "verbose",
        }
    },
    "loggers": {
        "road": {
            "handlers": ["file"],
            "level": "INFO",
            "propagate": True,
        }
    },
}

# SSL
CSRF_COOKIE_SECURE = True
SESSION_COOKIE_SECURE = True
SECURE_SSL_REDIRECT = True

# HSTS
# SECURE_HSTS_SECONDS = 60
# SECURE_HSTS_INCLUDE_SUBDOMAINS = True
# SECURE_HSTS_PRELOAD = True
