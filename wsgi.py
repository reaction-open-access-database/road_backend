"""
WSGI config for road project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.0/howto/deployment/wsgi/
"""

import os  # pylint: disable=unused-import

from django.core.wsgi import get_wsgi_application

application = get_wsgi_application()
