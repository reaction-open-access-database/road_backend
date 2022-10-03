from django.test import TestCase
from .models import Molecule, Reaction, ReactionComponent
from rdkit.Chem import AllChem
from django_rdkit.config import config
from .services import reaction_create, get_reactions_for_molecule
from django.contrib.auth.models import User


class MoleculeTest(TestCase):
    def setUp(self) -> None:
        self._user = User.objects.create_user('test')

    def test_simple_smiles_molecules(self):
        molecules = (
            'C',
            'CCl',
            'c1ccccc1',
        )

        for molecule in molecules:
            Molecule.objects.create(molecule=AllChem.MolFromSmiles(molecule), owner=self._user)

        for molecule in molecules:
            self.assertEqual(1, Molecule.objects.filter(molecule=molecule).count())

    def test_chiral_molecules(self):
        config.do_chiral_sss = True
        molecules = (
            # InChI
            AllChem.MolFromInchi('InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'),  # l-alanine
            AllChem.MolFromInchi('InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1'),  # r-alanine
            # SMILES
            AllChem.MolFromSmiles('C[C@@H](N)O'),  # (1S)-1-Aminoethanol
        )

        for molecule in molecules:
            Molecule.objects.create(molecule=molecule, owner=self._user)

        # If this fails, it's most likely because stereochemistry is not working
        for molecule in molecules:
            self.assertEqual(1, Molecule.objects.filter(molecule=molecule).count())

    # def test_molecular_formula(self):
    #     pass

    # def test_inchi_molecules(self):
    #     pass

    # def test_numbered_molecules(self):
    #     pass


class ReactionTest(TestCase):
    def setUp(self) -> None:
        self._user = User.objects.create_user('test')
    # def test_simple_reaction(self):
    #     reaction = AllChem.ReactionFromSmarts('CCl>[OH-]>CO')
    #     reaction_create(reaction)

    # def test_applying_reaction(self):
    #     pass

    def test_find_reactions(self):
        sample_reactions = (
            AllChem.ReactionFromSmarts('CCl>[OH-]>CO'),
            AllChem.ReactionFromSmarts('C.O=O>>C(=O)=O'),
        )

        for reaction in sample_reactions:
            reaction_create(reaction, self._user)

        self.assertEqual(2, Reaction.objects.count())
        self.assertEqual(1, get_reactions_for_molecule('CCl').count())
        self.assertEqual(1, get_reactions_for_molecule('CCl', ReactionComponent.ComponentType.REACTANT).count())
        self.assertEqual(0, get_reactions_for_molecule('CCl', ReactionComponent.ComponentType.PRODUCT).count())
        self.assertEqual(0, get_reactions_for_molecule('CCl', ReactionComponent.ComponentType.AGENT).count())
        self.assertEqual(1, get_reactions_for_molecule('CO').count())
        self.assertEqual(1, get_reactions_for_molecule('[OH-]').count())


class UserAccountTest(TestCase):
    def setUp(self) -> None:
        self._user = User.objects.create_user('test')
    # def test_create_user(self):
    #     pass
    pass
