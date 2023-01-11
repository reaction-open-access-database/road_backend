"""
The permissions for the ROAD REST API.
"""

from typing import Any

from rest_framework.permissions import SAFE_METHODS, BasePermission
from rest_framework.request import Request
from rest_framework.views import APIView

# Ensure that every class contains has_permission AND
# has_object_permission methods


class ReadOnly(BasePermission):
    """Allows access only to read-only methods."""

    def has_permission(self, request: Request, view: APIView) -> bool:
        return request.method in SAFE_METHODS

    def has_object_permission(self, request: Request, view: APIView, obj: Any) -> bool:
        return self.has_permission(request, view)
