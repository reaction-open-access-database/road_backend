from rest_framework.exceptions import APIException


class InvalidMolecule(APIException):
    status_code = 400
    default_detail = 'Invalid molecule'
    default_code = 'invalid_molecule'


class InvalidReaction(APIException):
    status_code = 400
    default_detail = 'Invalid reaction'
    default_code = 'invalid_reaction'
