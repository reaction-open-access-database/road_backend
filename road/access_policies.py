"""
The access policies to define the permissions for ROAD.
"""

from typing import Any

from rest_access_policy.access_policy import AccessPolicy
from rest_access_policy.access_view_set_mixin import AccessViewSetMixin

# Safe actions are: "list", "retrieve", "metadata"
# Unsafe actions are: "create", "update", "partial_update", "destroy"

SUPERUSER_ALLOW_ALL = {
    "action": "*",
    "principal": "admin",
    "effect": "allow",
}
OWNER_ALLOW_UPDATE = {
    "action": ["update", "partial_update", "destroy"],
    "principal": "authenticated",
    "effect": "allow",
    "condition": "is_owner",
}
AUTHENTICATED_ALLOW_CREATE = {
    "action": "create",
    "principal": "authenticated",
    "effect": "allow",
}
ANYONE_ALLOW_READ = {
    "action": ["list", "retrieve", "metadata"],
    "principal": "*",
    "effect": "allow",
}


class MoleculeAccessPolicy(AccessPolicy):
    """
    Access policy for the Molecule model.
    """

    statements = [
        SUPERUSER_ALLOW_ALL,
        OWNER_ALLOW_UPDATE,
        AUTHENTICATED_ALLOW_CREATE,
        ANYONE_ALLOW_READ,
    ]


class ReactionAccessPolicy(AccessPolicy):
    """
    Access policy for the Reaction model.
    """

    statements = [
        SUPERUSER_ALLOW_ALL,
    ]


class ReactionComponentAccessPolicy(AccessPolicy):
    """
    Access policy for the ReactionComponent model.
    """

    statements = [
        SUPERUSER_ALLOW_ALL,
    ]


class UserProfileAccessPolicy(AccessPolicy):
    """
    Access policy for the UserProfile model.
    """

    statements = [
        SUPERUSER_ALLOW_ALL,
        {
            "action": "retrieve",
            "principal": "authenticated",
            "effect": "allow",
            "condition": "is_owner",
        },
    ]


class OverrideAccessViewSetMixin(
    AccessViewSetMixin
):  # pylint: disable=too-few-public-methods
    """
    Mixin to set the permission classes for a ViewSet to the access policy.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)  # type: ignore

        self.permission_classes = [self.access_policy]
