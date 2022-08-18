from .models import Molecule
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
        fields = ['name', 'molecule']
