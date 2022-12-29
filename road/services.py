"""
Miscellaneous functions for creating and retrieving molecules and reactions.
"""

from typing import Optional

from django.contrib.auth.models import User  # pylint: disable=imported-auth-user
from django.db import transaction
from django.db.models import QuerySet
from rdkit.Chem.AllChem import Mol, ChemicalReaction

from .models import Molecule, Reaction, ReactionComponent


def get_reactions_for_molecule(
    molecule: Mol, component_type: Optional[ReactionComponent.ComponentType] = None
) -> QuerySet[Molecule]:
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
        return Reaction.objects.filter(components__molecule=molecule)  # type: ignore[no-any-return]

    return Reaction.objects.filter(  # type: ignore[no-any-return]
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
    molecule,
    component_type: ReactionComponent.ComponentType,
    owner: User,
) -> ReactionComponent:
    """
    Creates a reaction component for a reaction and molecule.
    Creates the molecule if it does not already exist.
    Returns the created reaction component.
    """
    molecule = molecule_get_or_create(molecule, owner)
    return ReactionComponent.objects.create(  # type: ignore[no-any-return]
        reaction=reaction,
        molecule=molecule,
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
    return molecule  # type: ignore[no-any-return]
