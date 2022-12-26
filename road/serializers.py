"""
The serializers for the ROAD REST API.
"""

import json

from django.contrib.auth.models import User
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rest_framework import serializers

from .models import Molecule, Reaction, ReactionComponent, UserProfile
from .exceptions import InvalidMolecule


class RDKitMoleculeJSONField(serializers.Field):
    """
    A field that serializes and deserializes RDKit molecules to and from JSON.
    """

    def to_representation(self, value):
        """Convert the RDKit molecule to JSON."""
        return json.loads(Chem.MolToJSON(value.molecule))

    def to_internal_value(self, data):
        """Convert the JSON to an RDKit molecule."""
        if data == "":
            return {"json": None}

        if isinstance(data, dict):
            json_data = json.dumps(data)
        else:
            json_data = data

        try:
            mols = Chem.JSONToMols(json_data)
        except RuntimeError:
            raise InvalidMolecule("Invalid JSON data")

        if len(mols) != 1:
            raise InvalidMolecule(
                f"Molecule field must contain exactly one molecule. "
                f"{len(mols)} molecules found."
            )

        return {"json": mols[0]}


class RDKitMoleculeSmilesField(serializers.Field):
    """
    A field that serializes and deserializes RDKit molecules to and from SMILES.
    """

    def to_representation(self, value):
        """Convert the RDKit molecule to SMILES."""
        return Chem.MolToSmiles(value.molecule)

    def to_internal_value(self, data):
        """Convert the SMILES to an RDKit molecule."""
        if data == "":
            return {"smiles": None}

        try:
            return {"smiles": Chem.MolFromSmiles(data)}
        except ValueError:
            raise InvalidMolecule("Invalid SMILES data")


class RDKitMoleculeInchiField(serializers.Field):
    """
    A field that serializes and deserializes RDKit molecules to and from InChI strings.
    """

    def to_representation(self, value):
        """Convert the RDKit molecule to an InChI string."""
        return Chem.MolToInchi(value.molecule)

    def to_internal_value(self, data):
        """Convert an InChI string to an RDKit molecule."""
        if data == "":
            return {"inchi": None}

        try:
            return {"inchi": Chem.MolFromInchi(data)}
        except ValueError:
            raise InvalidMolecule("Invalid InChI data")


class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer for the Molecule model.
    Displays the JSON, SMILES, InChI and SVG representations of the molecule.
    Also includes the name, molecular weight and molecular formula.
    """

    json = RDKitMoleculeJSONField(source="*", required=False)
    smiles = RDKitMoleculeSmilesField(source="*", required=False)
    inchi = RDKitMoleculeInchiField(source="*", required=False)
    svg = serializers.SerializerMethodField()
    mw = serializers.SerializerMethodField()
    formula = serializers.SerializerMethodField()
    id = serializers.ReadOnlyField()

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

    def validate(self, validated_data):
        validated_data.setdefault("json", None)
        validated_data.setdefault("smiles", None)
        validated_data.setdefault("inchi", None)

        representations = (
            validated_data.pop("json"),
            validated_data.pop("smiles"),
            validated_data.pop("inchi"),
        )

        provided_representations = [rep for rep in representations if rep]

        if len(provided_representations) == 1:
            validated_data["molecule"] = provided_representations[0]
        else:
            raise InvalidMolecule(
                "Exactly one of JSON, SMILES or InChI must be provided."
            )

        return validated_data

    def get_svg(self, obj):
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
        drawer.DrawMolecule(obj.molecule)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace("svg:", "")

    def get_mw(self, obj):
        return Chem.rdMolDescriptors.CalcExactMolWt(obj.molecule)

    def get_formula(self, obj):
        return Chem.rdMolDescriptors.CalcMolFormula(obj.molecule)


class ReactionSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer for the Reaction model.
    Displays the components of the reaction.
    """

    components = serializers.HyperlinkedRelatedField(
        many=True, read_only=True, view_name="reactioncomponent-detail"
    )

    class Meta:
        model = Reaction
        fields = ["url", "components"]


class ReactionComponentSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer for the ReactionComponent model.
    Displays the reaction, the molecule and the role of the molecule in the reaction.
    """

    class Meta:
        model = ReactionComponent
        fields = ["url", "reaction", "molecule", "component_type"]


class UserProfileSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer for the UserProfile model.
    Displays the user's username and email address.
    """

    username = serializers.ReadOnlyField(source="user.username")
    email = serializers.ReadOnlyField(source="user.email")

    class Meta:
        model = UserProfile
        fields = ["url", "username", "email"]


class UserSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer for the User model.
    Displays the user's profile.
    """

    profile = serializers.HyperlinkedRelatedField(
        read_only=True, view_name="userprofile-detail"
    )

    class Meta:
        model = User
        fields = ["profile"]
