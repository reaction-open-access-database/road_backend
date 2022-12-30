"""
The serializers for the ROAD REST API.
"""

import json
from typing import Any, Dict, Optional

from django.contrib.auth.models import User  # pylint: disable=imported-auth-user
from rdkit import Chem
from rdkit.Chem.AllChem import Mol
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcMolFormula
from rest_framework.serializers import (
    Field,
    HyperlinkedModelSerializer,
    HyperlinkedRelatedField,
    ReadOnlyField,
    SerializerMethodField,
)

from .exceptions import InvalidMolecule
from .models import Molecule, Reaction, ReactionComponent, UserProfile


class RDKitMoleculeJSONField(Field):  # type: ignore
    """
    A field that serializes and deserializes RDKit molecules to and from JSON.
    """

    def to_representation(self, value: Molecule) -> Any:
        """Convert the Molecule to JSON."""
        return json.loads(Chem.MolToJSON(value.molecule))

    def to_internal_value(self, data: str) -> Dict[str, Optional[Mol]]:
        """Convert the JSON to an RDKit molecule."""
        if data == "":
            return {"json": None}

        if isinstance(data, dict):
            json_data = json.dumps(data)
        else:
            json_data = data

        try:
            mols = Chem.JSONToMols(json_data)
        except RuntimeError as exc:
            raise InvalidMolecule("Invalid JSON data") from exc

        if len(mols) != 1:
            raise InvalidMolecule(
                f"Molecule field must contain exactly one molecule. "
                f"{len(mols)} molecules found."
            )

        return {"json": mols[0]}


class RDKitMoleculeSmilesField(Field):  # type: ignore
    """
    A field that serializes and deserializes RDKit molecules to and from SMILES.
    """

    def to_representation(self, value: Molecule) -> str:
        """Convert the RDKit molecule to SMILES."""
        return Chem.MolToSmiles(value.molecule)

    def to_internal_value(self, data: str) -> Dict[str, Optional[Mol]]:
        """Convert the SMILES to an RDKit molecule."""
        if data == "":
            return {"smiles": None}

        try:
            return {"smiles": Chem.MolFromSmiles(data)}
        except ValueError as exc:
            raise InvalidMolecule("Invalid SMILES data") from exc


class RDKitMoleculeInchiField(Field):  # type: ignore
    """
    A field that serializes and deserializes RDKit molecules to and from InChI strings.
    """

    def to_representation(self, value: Molecule) -> str:
        """Convert the RDKit molecule to an InChI string."""
        return Chem.MolToInchi(value.molecule)

    def to_internal_value(self, data: str) -> Dict[str, Optional[Mol]]:
        """Convert an InChI string to an RDKit molecule."""
        if data == "":
            return {"inchi": None}

        try:
            return {"inchi": Chem.MolFromInchi(data)}
        except ValueError as exc:
            raise InvalidMolecule("Invalid InChI data") from exc


class MoleculeSerializer(HyperlinkedModelSerializer):
    """
    Serializer for the Molecule model.
    Displays the JSON, SMILES, InChI and SVG representations of the molecule.
    Also includes the name, molecular weight and molecular formula.
    """

    json = RDKitMoleculeJSONField(source="*", required=False)
    smiles = RDKitMoleculeSmilesField(source="*", required=False)
    inchi = RDKitMoleculeInchiField(source="*", required=False)
    svg = SerializerMethodField()
    mw = SerializerMethodField()
    formula = SerializerMethodField()
    id = ReadOnlyField()

    class Meta:
        model = Molecule
        fields = [
            "url",
            "name",
            "json",
            "smiles",
            "inchi",
            "svg",
            "mw",
            "formula",
            "id",
        ]

    def validate(self, attrs: Dict[str, Any]) -> Dict[str, Any]:
        attrs.setdefault("json", None)
        attrs.setdefault("smiles", None)
        attrs.setdefault("inchi", None)

        representations = (
            attrs.pop("json"),
            attrs.pop("smiles"),
            attrs.pop("inchi"),
        )

        provided_representations = [rep for rep in representations if rep]

        if len(provided_representations) == 1:
            attrs["molecule"] = provided_representations[0]
        else:
            raise InvalidMolecule(
                "Exactly one of JSON, SMILES or InChI must be provided."
            )

        return attrs

    def get_svg(self, obj: Molecule) -> str:
        """Return the text of an SVG image of the molecule."""
        drawer = MolDraw2DSVG(200, 200)
        drawer.DrawMolecule(obj.molecule)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace("svg:", "")

    def get_mw(self, obj: Molecule) -> float:
        """Return the molecular weight of the molecule."""
        return CalcExactMolWt(obj.molecule)

    def get_formula(self, obj: Molecule) -> str:
        """Return the molecular formula of the molecule."""
        return CalcMolFormula(obj.molecule)


class ReactionSerializer(HyperlinkedModelSerializer):
    """
    Serializer for the Reaction model.
    Displays the components of the reaction.
    """

    components: HyperlinkedRelatedField = HyperlinkedRelatedField(  # type: ignore
        many=True, read_only=True, view_name="reactioncomponent-detail"
    )

    class Meta:
        model = Reaction
        fields = ["url", "components"]


class ReactionComponentSerializer(HyperlinkedModelSerializer):
    """
    Serializer for the ReactionComponent model.
    Displays the reaction, the molecule and the role of the molecule in the reaction.
    """

    class Meta:
        model = ReactionComponent
        fields = ["url", "reaction", "molecule", "component_type"]


class UserProfileSerializer(HyperlinkedModelSerializer):
    """
    Serializer for the UserProfile model.
    Displays the user's username and email address.
    """

    username = ReadOnlyField(source="owner.username")
    email = ReadOnlyField(source="owner.email")

    class Meta:
        model = UserProfile
        fields = ["url", "username", "email"]


class UserSerializer(HyperlinkedModelSerializer):
    """
    Serializer for the User model.
    Displays the user's profile.
    """

    profile: HyperlinkedRelatedField = HyperlinkedRelatedField(  # type: ignore
        read_only=True, view_name="userprofile-detail"
    )

    class Meta:
        model = User
        fields = ["profile"]
