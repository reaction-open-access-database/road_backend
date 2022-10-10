from .models import Molecule, Reaction, ReactionComponent
from rdkit.Chem import AllChem
from django_rdkit.config import config
from .services import reaction_create, get_reactions_for_molecule
from django.contrib.auth.models import User
from rest_framework.test import APITestCase
from django.urls import reverse
from rest_framework import status


class MoleculeTest(APITestCase):
    def setUp(self) -> None:
        self._user = User.objects.create_user('test')
        self.client.force_authenticate(user=self._user)

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

    def test_create_molecule_smiles(self):
        smiles = 'O=O'
        name = 'oxygen'
        response = self.client.post(reverse('molecule-list'), {'name': name, 'smiles': smiles}, format='json')

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(AllChem.MolToSmiles(Molecule.objects.get().molecule), smiles)
        self.assertEqual(Molecule.objects.get().name, 'oxygen')
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_inchi(self):
        inchi = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'
        name = 'l-alanine'
        response = self.client.post(reverse('molecule-list'), {'name': name, 'inchi': inchi}, format='json')

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(AllChem.MolToInchi(Molecule.objects.get().molecule), inchi)
        self.assertEqual(Molecule.objects.get().name, 'l-alanine')
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_json(self):
        molecule = AllChem.MolFromSmiles('C=C-C')
        json = AllChem.MolToJSON(molecule)
        name = 'propene'
        response = self.client.post(reverse('molecule-list'), {'name': name, 'json': json}, format='json')

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(AllChem.MolToInchi(Molecule.objects.get().molecule), AllChem.MolToInchi(molecule))
        self.assertEqual(Molecule.objects.get().name, 'propene')
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_nothing(self):
        name = 'nothing'
        response = self.client.post(reverse('molecule-list'), {'name': name}, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(Molecule.objects.count(), 0)

    # def test_molecular_formula(self):
    #     pass

    # def test_inchi_molecules(self):
    #     pass

    # def test_numbered_molecules(self):
    #     pass


class ReactionTest(APITestCase):
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


class UserAccountTest(APITestCase):
    def setUp(self) -> None:
        self._user = User.objects.create_user('test')
    # def test_create_user(self):
    #     pass
    pass
