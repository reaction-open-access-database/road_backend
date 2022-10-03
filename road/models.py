from django_rdkit import models
from rdkit import Chem
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver


class Molecule(models.Model):
    name = models.CharField(max_length=256)
    molecule = models.MolField()
    owner = models.ForeignKey(User, on_delete=models.RESTRICT,
                              related_name='molecules')
    molecular_formula = models.CharField(max_length=256)

    def save(self, *args, **kwargs):
        self.molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(self.molecule)
        super().save(*args, **kwargs)

    def get_inchi(self):
        return Chem.MolToInchi(self.molecule)

    def get_smiles(self):
        return Chem.MolToSmiles(self.molecule)


class Reaction(models.Model):
    owner = models.ForeignKey(User, on_delete=models.RESTRICT,
                              related_name='reactions')


class ReactionComponent(models.Model):
    class ComponentType(models.TextChoices):
        REACTANT = 'reactant'
        PRODUCT = 'product'
        AGENT = 'agent'

    component_type = models.CharField(
        choices=ComponentType.choices,
        max_length=10,
    )
    reaction = models.ForeignKey(Reaction, on_delete=models.CASCADE,
                                 related_name='components')
    molecule = models.ForeignKey(Molecule, on_delete=models.RESTRICT,
                                 related_name='components')
    count_numerator = models.IntegerField(default=1)
    count_denominator = models.IntegerField(default=1)
    owner = models.ForeignKey(User, on_delete=models.RESTRICT,
                              related_name='components')

    def get_component_type(self):
        return self.ComponentType[self.component_type]


class ReactionSource(models.Model):
    """
    Contains information about where the reaction came from.
    This could be a journal article, book, website or other source.
    """
    reaction = models.ForeignKey(Reaction, on_delete=models.CASCADE,
                                 related_name='reaction_sources')
    owner = models.ForeignKey(User, on_delete=models.RESTRICT,
                              related_name='reaction_sources')


class UserProfile(models.Model):
    owner = models.OneToOneField(User, on_delete=models.CASCADE,
                                 related_name='profile')

    def __str__(self):
        return self.owner.username


# Automatically create UserProfile when a User is created
@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        UserProfile.objects.create(owner=instance)
