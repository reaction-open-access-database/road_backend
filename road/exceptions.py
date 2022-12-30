"""
The exceptions used in the ROAD REST API.
"""

from rest_framework.exceptions import APIException


class InvalidMolecule(APIException):
    """The molecule is invalid."""

    status_code = 400
    default_detail = "Invalid molecule"
    default_code = "invalid_molecule"


class InvalidReaction(APIException):
    """The reaction is invalid."""

    status_code = 400
    default_detail = "Invalid reaction"
    default_code = "invalid_reaction"


class InvalidQuery(APIException):
    """The query is invalid."""

    status_code = 400
    default_detail = "Invalid query"
    default_code = "invalid_query"


class ParameterNotProvided(APIException):
    """A required parameter was not provided."""

    status_code = 400
    default_detail = "Argument not provided"
    default_code = "argument_not_provided"
