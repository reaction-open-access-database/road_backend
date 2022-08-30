from .models import Molecule, UserProfile
from .exceptions import InvalidMolecule
from rest_framework import serializers
from rdkit import Chem
import json


class RDKitMoleculeField(serializers.Field):
    def to_representation(self, value):
        return json.loads(Chem.MolToJSON(value))

    def to_internal_value(self, data):
        json_data = json.dumps(data)

        try:
            mols = Chem.JSONToMols(json_data)
        except RuntimeError:
            raise InvalidMolecule('Invalid JSON data')

        if len(mols) != 1:
            raise InvalidMolecule(
                f'Molecule field must contain exactly one molecule. '
                f'{len(mols)} molecules found.'
            )

        return mols[0]


class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    molecule = RDKitMoleculeField()

    class Meta:
        model = Molecule
        fields = ['url', 'name', 'molecule']


class UserProfileSerializer(serializers.HyperlinkedModelSerializer):
    username = serializers.ReadOnlyField(source='user.username')
    email = serializers.ReadOnlyField(source='user.email')

    class Meta:
        model = UserProfile
        fields = ['url', 'username', 'email']
