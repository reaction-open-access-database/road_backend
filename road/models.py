from django_rdkit import models
from rdkit import Chem


class Molecule(models.Model):
    name = models.CharField(max_length=256)
    molecule = models.MolField()

    def get_inchi(self):
        return Chem.MolToInchi(self.molecule)

    def get_smiles(self):
        return Chem.MolToSmiles(self.molecule)


class Reaction(models.Model):
    pass


class ReactionComponent(models.Model):
    class ComponentType(models.TextChoices):
        REACTANT = 'reactant'
        PRODUCT = 'product'
        AGENT = 'agent'

    component_type = models.CharField(
        choices=ComponentType.choices,
        max_length=10,
    )
    reaction = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    count = models.IntegerField(default=1)  # TODO: Allow fractions

    def get_component_type(self):
        return self.ComponentType[self.component_type]
