"""
The permissions for the ROAD REST API.
"""

from rest_framework.permissions import *  # pylint: disable=wildcard-import, unused-wildcard-import


# Ensure that every class contains has_permission AND
# has_object_permission methods


class IsOwner(BasePermission):
    """
    Allows the owner of an object to view and edit it,
    but other authenticated users to view it only.
    """

    def has_permission(self, request, view):
        return request.user and request.user.is_authenticated

    def has_object_permission(self, request, view, obj):
        return (
            request.user and request.user.is_authenticated and request.user == obj.owner
        )


class IsSuperUser(BasePermission):
    """Allows access only to superusers."""

    def has_permission(self, request, view):
        return (
            request.user and request.user.is_authenticated and request.user.is_superuser
        )

    def has_object_permission(self, request, view, obj):
        return self.has_permission(request, view)


class ReadOnly(BasePermission):
    """Allows access only to read-only methods."""

    def has_permission(self, request, view):
        return request.method in SAFE_METHODS

    def has_object_permission(self, request, view, obj):
        return self.has_permission(request, view)
