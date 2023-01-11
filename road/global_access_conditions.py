"""
The global access conditions used for the permission system.
"""

from typing import Any

from rest_framework.request import Request


def is_owner(
    request: Request, view: Any, action: str  # pylint: disable=unused-argument
) -> bool:
    """Returns whether the request's user is the owner of the object."""
    return bool(request.user == view.get_object().owner)
