import inspect

from typing import Type

from rest_access_policy import AccessPolicy
from rest_framework.response import Response

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
    statements = [
        SUPERUSER_ALLOW_ALL,
        OWNER_ALLOW_UPDATE,
        AUTHENTICATED_ALLOW_CREATE,
        ANYONE_ALLOW_READ,
    ]


class ReactionAccessPolicy(AccessPolicy):
    statements = [
        SUPERUSER_ALLOW_ALL,
    ]


class ReactionComponentAccessPolicy(AccessPolicy):
    statements = [
        SUPERUSER_ALLOW_ALL,
    ]


class UserProfileAccessPolicy(AccessPolicy):
    statements = [
        SUPERUSER_ALLOW_ALL,
        {
            "action": "retrieve",
            "principal": "authenticated",
            "effect": "allow",
            "condition": "is_owner",
        },
    ]


class OverrideAccessViewSetMixin(object):
    access_policy: Type[AccessPolicy]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        access_policy = getattr(self, "access_policy", None)

        if not inspect.isclass(access_policy) or not issubclass(
            access_policy, AccessPolicy
        ):
            raise Exception(
                """
                    When mixing AccessViewSetMixin into your view set, you must assign an AccessPolicy 
                    to the access_policy class attribute.
                """
            )

        self.permission_classes = [access_policy]

    def finalize_response(self, request, response, *args, **kwargs) -> Response:
        response = super().finalize_response(request, response, *args, **kwargs)
        return response
