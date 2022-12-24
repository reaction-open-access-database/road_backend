from typing import Optional

from django.contrib.auth.models import User
from django.db import transaction
from rdkit.Chem import AllChem

from .models import Molecule, Reaction, ReactionComponent


def get_reactions_for_molecule(
        molecule,
        component_type: Optional[ReactionComponent.ComponentType] = None):
    try:
        molecule = Molecule.objects.get(molecule=molecule)
    except Molecule.DoesNotExist:
        return []

    if component_type is None:
        return Reaction.objects.filter(components__molecule=molecule)
    else:
        return Reaction.objects.filter(
            components__molecule=molecule,
            components__component_type=component_type
        )


@transaction.atomic
def reaction_create(rdkit_reaction: AllChem.ChemicalReaction, owner: User) -> Reaction:
    reaction = Reaction.objects.create(owner=owner)

    # Create the reaction components
    for reactant in rdkit_reaction.GetReactants():
        reaction_component_create(reaction, reactant,
                                  ReactionComponent.ComponentType.REACTANT, owner)

    for product in rdkit_reaction.GetProducts():
        reaction_component_create(reaction, product,
                                  ReactionComponent.ComponentType.PRODUCT, owner)

    for agent in rdkit_reaction.GetAgents():
        reaction_component_create(reaction, agent,
                                  ReactionComponent.ComponentType.AGENT, owner)

    return reaction


def reaction_component_create(reaction: Reaction, molecule,
                              component_type: ReactionComponent.ComponentType,
                              owner: User) -> ReactionComponent:
    molecule = molecule_get_or_create(molecule, owner)
    return ReactionComponent.objects.create(
        reaction=reaction,
        molecule=molecule,
        component_type=component_type,
        owner=owner,
    )


def molecule_get_or_create(rdkit_molecule, owner: User) -> Molecule:
    try:
        molecule = Molecule.objects.get(molecule=rdkit_molecule)
    except Molecule.DoesNotExist:
        molecule = Molecule.objects.create(molecule=rdkit_molecule, owner=owner)
    return molecule
