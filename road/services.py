"""
Miscellaneous functions for creating and retrieving molecules and reactions.
"""

from typing import Any, Dict, Optional

from django.contrib.auth.models import User  # pylint: disable=imported-auth-user
from django.db import IntegrityError, transaction
from django.db.models import QuerySet
from rdkit.Chem.AllChem import ChemicalReaction, Mol
from rest_framework.response import Response
from rest_framework.status import HTTP_400_BAD_REQUEST
from rest_framework.views import exception_handler

from .models import Molecule, Reaction, ReactionComponent


def get_reactions_for_molecule(
    molecule: Mol, component_type: Optional[ReactionComponent.ComponentType] = None
) -> QuerySet[Reaction]:
    """
    Returns the reactions associated with a particular molecule.
    Will raise a Molecule.DoesNotExist exception if the molecule does not exist in the database.

    :param molecule: The molecule to check
    :param component_type: The reaction component type the molecule should be
    (reactant, product or agent), or None (default) to permit any component type
    :return: A list of the associated reactions
    """
    molecule = Molecule.objects.get(molecule=molecule)

    if component_type is None:
        return Reaction.objects.filter(components__molecule=molecule)

    return Reaction.objects.filter(
        components__molecule=molecule, components__component_type=component_type
    )


@transaction.atomic
def reaction_create(rdkit_reaction: ChemicalReaction, owner: User) -> Reaction:
    """
    Creates a reaction from an RDKit reaction, with the specified owner.
    Returns the created reaction.
    """
    reaction: Reaction = Reaction.objects.create(owner=owner)

    # Create the reaction components
    for reactant in rdkit_reaction.GetReactants():
        reaction_component_create(
            reaction, reactant, ReactionComponent.ComponentType.REACTANT, owner
        )

    for product in rdkit_reaction.GetProducts():
        reaction_component_create(
            reaction, product, ReactionComponent.ComponentType.PRODUCT, owner
        )

    for agent in rdkit_reaction.GetAgents():
        reaction_component_create(
            reaction, agent, ReactionComponent.ComponentType.AGENT, owner
        )

    return reaction


def reaction_component_create(
    reaction: Reaction,
    rdkit_molecule: Mol,
    component_type: ReactionComponent.ComponentType,
    owner: User,
) -> ReactionComponent:
    """
    Creates a reaction component for a reaction and molecule.
    Creates the molecule if it does not already exist.
    Returns the created reaction component.
    """
    rdkit_molecule = molecule_get_or_create(rdkit_molecule, owner)
    return ReactionComponent.objects.create(
        reaction=reaction,
        molecule=rdkit_molecule,
        component_type=component_type,
        owner=owner,
    )


def molecule_get_or_create(rdkit_molecule: Mol, owner: User) -> Molecule:
    """
    Creates a molecule if it does not already exist.
    Returns the molecule model.
    """
    try:
        molecule = Molecule.objects.get(molecule=rdkit_molecule)
    except Molecule.DoesNotExist:
        molecule = Molecule.objects.create(molecule=rdkit_molecule, owner=owner)
    return molecule


def custom_exception_handler(
    exc: Exception, context: Dict[str, Any]
) -> Optional[Response]:
    """A custom exception handler to return a 400 status code for database IntegrityErrors."""
    # From https://stackoverflow.com/a/50776557
    # Call REST framework's default exception handler first to get the standard error response.
    response = exception_handler(exc, context)
    print(context)

    # if there is an IntegrityError and the error response hasn't already been generated
    if isinstance(exc, IntegrityError) and not response:
        response = Response(
            {
                "message": "There is a conflict between the data you are trying to save "
                "and the current data."
            },
            status=HTTP_400_BAD_REQUEST,
        )

    return response
