"""
ASGI config for road project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.0/howto/deployment/asgi/
"""

import os  # pylint: disable=unused-import

from django.core.asgi import get_asgi_application

application = get_asgi_application()
