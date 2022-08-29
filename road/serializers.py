from .models import Molecule, UserProfile
from rest_framework import serializers
from rdkit import Chem
import json


class RDKitMoleculeField(serializers.Field):
    def to_representation(self, value):
        return json.loads(Chem.MolToJSON(value))

    def to_internal_value(self, data):
        mols = Chem.JSONToMols(json.dumps(data))
        assert len(mols) == 1
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
