from typing import Optional
from rdkit.Chem import AllChem
from django.db import transaction
from .models import Molecule, Reaction, ReactionComponent


def get_reactions_for_molecule(
        molecule,
        component_type: Optional[ReactionComponent.ComponentType]):
    try:
        molecule = Molecule.objects.get(molecule=molecule)
    except Molecule.DoesNotExist:
        return []

    if component_type is None:
        return Reaction.objects.filter(reactioncomponent__molecule=molecule)
    else:
        return Reaction.objects.filter(
            reactioncomponent__molecule=molecule,
            reactioncomponent__component_type=component_type
        )


@transaction.atomic
def reaction_create(rdkit_reaction: AllChem.ChemicalReaction) -> Reaction:
    reaction = Reaction.objects.create()

    # Create the reaction components
    for reactant in rdkit_reaction.GetReactants():
        reaction_component_create(reaction, reactant,
                                  ReactionComponent.ComponentType.REACTANT)

    for product in rdkit_reaction.GetProducts():
        reaction_component_create(reaction, product,
                                  ReactionComponent.ComponentType.PRODUCT)

    for agent in rdkit_reaction.GetAgents():
        reaction_component_create(reaction, agent,
                                  ReactionComponent.ComponentType.AGENT)

    # TODO: Connect the reaction component atoms

    return reaction


def reaction_component_create(reaction: Reaction, molecule, component_type: ReactionComponent.ComponentType) -> ReactionComponent:
    molecule = molecule_get_or_create(molecule)
    return ReactionComponent.objects.create(
        reaction=reaction,
        molecule=molecule,
        component_type=component_type,
    )


def molecule_get_or_create(rdkit_molecule) -> Molecule:
    try:
        molecule = Molecule.objects.get(molecule=rdkit_molecule)
    except Molecule.DoesNotExist:
        molecule = Molecule.objects.create(molecule=rdkit_molecule)
    return molecule
