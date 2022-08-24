from rest_framework.permissions import *


class IsOwnerOrReadOnly(BasePermission):
    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in SAFE_METHODS:
            return True

        # Write permissions are only allowed to the owner of the snippet.
        try:  # TODO: Change this to disallow other users from editing, even if there is no owner
            object_owner = obj.owner
        except AttributeError:
            return False

        return object_owner == request.user
