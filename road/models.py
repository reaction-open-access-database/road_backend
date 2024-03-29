"""
The models for ROAD.
"""

import json
from typing import Any, Iterable, Optional

from django.contrib.auth.models import User  # pylint: disable=imported-auth-user
from django.db.models.signals import post_save
from django.dispatch import receiver
from django_rdkit import models
from rdkit.Chem import JSONToMols, Mol, MolToInchi, MolToJSON, MolToSmiles, SanitizeMol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


class SerializableMolField(models.MolField):
    """
    A field that stores RDKit molecules in the database.
    Additionally, it provides the ability to serialize and deserialize molecules to/from JSON.
    """

    def value_to_string(self, obj: models.Model) -> str:
        mol = super().value_from_object(obj)  # type: ignore
        return MolToJSON(mol)

    def value_from_object(self, obj: models.Model) -> Any:
        return json.loads(self.value_to_string(obj))

    def to_python(self, value: Any) -> Mol:
        try:
            if isinstance(value, str):
                json_data = json.loads(value)
            else:
                json_data = value
            return JSONToMols(json.dumps(json_data))[0]
        except RuntimeError:
            return super().to_python(value)  # type: ignore


class Molecule(models.Model):
    """
    A molecule used in one or more reactions.
    The molecule is stored as an RDKit molecule.
    The molecular formula is derived from the RDKit molecule, so should not be written directly.
    """

    name = models.CharField(max_length=256)
    molecule = SerializableMolField(unique=True)  # type: ignore
    owner = models.ForeignKey(User, on_delete=models.RESTRICT, related_name="molecules")
    molecular_formula = models.CharField(max_length=256)

    def save(
        self,
        force_insert: bool = False,
        force_update: bool = False,
        using: Optional[str] = None,
        update_fields: Optional[Iterable[str]] = None,
    ) -> None:
        """Save the molecule."""
        SanitizeMol(self.molecule)
        self.molecular_formula = CalcMolFormula(self.molecule)
        super().save(force_insert, force_update, using, update_fields)

    def get_inchi(self) -> str:
        """Return the InChI representation of the molecule."""
        return MolToInchi(self.molecule)

    def get_smiles(self) -> str:
        """Return the SMILES representation of the molecule."""
        return MolToSmiles(self.molecule)


class Reaction(models.Model):
    """
    A reaction from one or more reactants to one or more products, possibly with some agents.
    """

    owner = models.ForeignKey(User, on_delete=models.RESTRICT, related_name="reactions")


class ReactionComponent(models.Model):
    """
    An individual component of a reaction.
    This could be a reactant, product, or agent (e.g. catalyst, solvent).
    Also allows for fractional amounts using the count_numerator and count_denominator fields.
    """

    class ComponentType(models.TextChoices):
        """
        The role of the component in the reaction.
        One of a reactant, product, or agent.
        """

        REACTANT = "reactant"
        PRODUCT = "product"
        AGENT = "agent"

    component_type = models.CharField(
        choices=ComponentType.choices,
        max_length=10,
    )
    reaction = models.ForeignKey(
        Reaction, on_delete=models.CASCADE, related_name="components"
    )
    molecule = models.ForeignKey(
        Molecule, on_delete=models.RESTRICT, related_name="components"
    )
    count_numerator = models.IntegerField(default=1)
    count_denominator = models.IntegerField(default=1)
    owner = models.ForeignKey(
        User, on_delete=models.RESTRICT, related_name="components"
    )

    def get_component_type(self) -> str:
        """Return the component type (reactant, product or agent) as a string."""
        return self.ComponentType[self.component_type]


class ReactionSource(models.Model):
    """
    Contains information about where the reaction came from.
    This could be a journal article, book, website or other source.
    """

    reaction = models.ForeignKey(
        Reaction, on_delete=models.CASCADE, related_name="reaction_sources"
    )
    owner = models.ForeignKey(
        User, on_delete=models.RESTRICT, related_name="reaction_sources"
    )


class UserProfile(models.Model):
    """
    Stores extra information about a user.
    """

    owner = models.OneToOneField(User, on_delete=models.CASCADE, related_name="profile")

    def __str__(self) -> str:
        return self.owner.username


# Automatically create UserProfile when a User is created
@receiver(post_save, sender=User)
def create_user_profile(
    instance: User,
    created: bool,
    **kwargs: Any,
) -> None:
    """Create a UserProfile when a User is created."""
    if created:
        UserProfile.objects.create(owner=instance)
