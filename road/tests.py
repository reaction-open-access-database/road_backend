"""
Tests ROAD.

MoleculeTest:
    Tests that the Molecule model is correct.
ReactionTest:
    Tests that the Reaction model is correct.
"""

from django.contrib.auth.models import (  # pylint: disable=imported-auth-user
    User,
    AnonymousUser,
)
from django.urls import reverse
from django_rdkit.config import config
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles
from rest_framework import status
from rest_framework.test import APITestCase

from .models import Molecule, Reaction, ReactionComponent
from .services import get_reactions_for_molecule, reaction_create


class MoleculeTest(APITestCase):
    """Tests that the Molecule model is correct."""

    def setUp(self) -> None:
        """Set up the dummy user for the tests."""
        self._user = User.objects.create_user("test")
        self.client.force_authenticate(user=self._user)

    def test_simple_smiles_molecules(self) -> None:
        """Test that simple SMILES molecules are correctly parsed and stored."""
        molecules = (
            "C",
            "CCl",
            "c1ccccc1",
        )

        for molecule in molecules:
            Molecule.objects.create(
                molecule=AllChem.MolFromSmiles(molecule), owner=self._user
            )

        for molecule in molecules:
            self.assertEqual(1, Molecule.objects.filter(molecule=molecule).count())

    def test_chiral_molecules(self) -> None:
        """Test that chiral molecules are correctly parsed and stored."""
        config.do_chiral_sss = True
        molecules = (
            # InChI
            AllChem.MolFromInchi(
                "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"
            ),  # l-alanine
            AllChem.MolFromInchi(
                "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1"
            ),  # r-alanine
            # SMILES
            AllChem.MolFromSmiles("C[C@@H](N)O"),  # (1S)-1-Aminoethanol
        )

        for molecule in molecules:
            Molecule.objects.create(molecule=molecule, owner=self._user)

        # If this fails, it's most likely because stereochemistry is not working
        for molecule in molecules:
            self.assertEqual(1, Molecule.objects.filter(molecule=molecule).count())

    def test_create_molecule_smiles(self) -> None:
        """Test that a molecule can be created from a SMILES string through the REST API."""
        smiles = "O=O"
        name = "oxygen"
        response = self.client.post(
            reverse("molecule-list"), {"name": name, "smiles": smiles}, format="json"
        )

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(AllChem.MolToSmiles(Molecule.objects.get().molecule), smiles)
        self.assertEqual(Molecule.objects.get().name, "oxygen")
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_inchi(self) -> None:
        """Test that a molecule can be created from an InChI string through the REST API."""
        inchi = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"
        name = "l-alanine"
        response = self.client.post(
            reverse("molecule-list"), {"name": name, "inchi": inchi}, format="json"
        )

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(AllChem.MolToInchi(Molecule.objects.get().molecule), inchi)
        self.assertEqual(Molecule.objects.get().name, "l-alanine")
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_json(self) -> None:
        """Test that a molecule can be created from a JSON string through the REST API."""
        molecule = AllChem.MolFromSmiles("C=C-C")
        json = AllChem.MolToJSON(molecule)
        name = "propene"
        response = self.client.post(
            reverse("molecule-list"), {"name": name, "json": json}, format="json"
        )

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Molecule.objects.count(), 1)
        self.assertEqual(
            AllChem.MolToInchi(Molecule.objects.get().molecule),
            AllChem.MolToInchi(molecule),
        )
        self.assertEqual(Molecule.objects.get().name, "propene")
        self.assertEqual(Molecule.objects.get().owner, self._user)

    def test_create_molecule_nothing(self) -> None:
        """Test that a molecule cannot be created without any data through the REST API."""
        name = "nothing"
        response = self.client.post(
            reverse("molecule-list"), {"name": name}, format="json"
        )

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(Molecule.objects.count(), 0)

    def test_serializable_mol_field(self) -> None:
        """
        Test that the SerializableMolField class can correctly serialize
        and deserialize molecules.
        """
        mol = AllChem.MolFromSmiles("C=C-C")
        name = "propene"

        molecule = Molecule.objects.create(molecule=mol, name=name, owner=self._user)

        db_json_str = Molecule.molecule.field.value_to_string(molecule)  # type: ignore

        db_molecule = AllChem.JSONToMols(db_json_str)[0]

        self.assertTrue(mol.HasSubstructMatch(db_molecule))
        self.assertTrue(db_molecule.HasSubstructMatch(mol))

    # def test_molecular_formula(self):
    #     pass

    # def test_inchi_molecules(self):
    #     pass

    # def test_numbered_molecules(self):
    #     pass


class ReactionTest(APITestCase):
    """Tests that the Reaction model is correct."""

    def setUp(self) -> None:
        """Set up the dummy user for these tests."""
        self._user = User.objects.create_user("test")

    # def test_simple_reaction(self):
    #     reaction = AllChem.ReactionFromSmarts('CCl>[OH-]>CO')
    #     reaction_create(reaction)

    # def test_applying_reaction(self):
    #     pass

    def test_find_reactions(self) -> None:
        """Test finding reactions that contain a molecule."""
        sample_reactions = (
            AllChem.ReactionFromSmarts("CCl>[OH-]>CO"),
            AllChem.ReactionFromSmarts("C.O=O>>C(=O)=O"),
        )

        for reaction in sample_reactions:
            reaction_create(reaction, self._user)

        self.assertEqual(2, Reaction.objects.count())
        self.assertEqual(1, get_reactions_for_molecule("CCl").count())
        self.assertEqual(
            1,
            get_reactions_for_molecule(
                "CCl", ReactionComponent.ComponentType.REACTANT
            ).count(),
        )
        self.assertEqual(
            0,
            get_reactions_for_molecule(
                "CCl", ReactionComponent.ComponentType.PRODUCT
            ).count(),
        )
        self.assertEqual(
            0,
            get_reactions_for_molecule(
                "CCl", ReactionComponent.ComponentType.AGENT
            ).count(),
        )
        self.assertEqual(1, get_reactions_for_molecule("CO").count())
        self.assertEqual(1, get_reactions_for_molecule("[OH-]").count())


class PermissionTest(APITestCase):
    """Tests that the permissions are correct."""

    def setUp(self) -> None:
        """Set up the dummy anonymous, authenticated and admin users for these tests."""
        self._normal_user = User.objects.create_user("normal")
        self._admin_user = User.objects.create_superuser("admin")
        self._anonymous_user = AnonymousUser()

        self._normal_molecule = Molecule.objects.create(
            name="methane", molecule=MolFromSmiles("C"), owner=self._normal_user
        )
        self._admin_molecule = Molecule.objects.create(
            name="ethane", molecule=MolFromSmiles("CC"), owner=self._admin_user
        )

    def test_molecule_anonymous_permissions(self) -> None:
        """Test that anonymous users can only read molecules."""
        self.client.force_authenticate(user=self._anonymous_user)

        with self.subTest("Anonymous users are not allowed to create molecules"):
            response = self.client.post(
                reverse("molecule-list"),
                {"name": "methane", "smiles": "C"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)

        with self.subTest("Anonymous users can list molecules"):
            response = self.client.get(reverse("molecule-list"))
            self.assertEqual(response.status_code, status.HTTP_200_OK)

        with self.subTest("Anonymous users can retrieve molecules"):
            response = self.client.get(
                reverse("molecule-detail", args=[self._normal_molecule.id])
            )
            self.assertEqual(response.status_code, status.HTTP_200_OK)

        with self.subTest("Anonymous users cannot update molecules"):
            response = self.client.put(
                reverse("molecule-detail", args=[self._normal_molecule.id]),
                {"name": "test", "smiles": "C"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
            self.assertEqual(self._normal_molecule.name, "methane")

    def test_molecule_authenticated_permissions(self) -> None:
        """Test that authenticated users can only read and update their own molecules."""
        self.client.force_authenticate(user=self._normal_user)

        with self.subTest("Authenticated users are allowed to create molecules"):
            response = self.client.post(
                reverse("molecule-list"),
                {"name": "methane", "smiles": "C"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_201_CREATED)

        with self.subTest(
            "Authenticated users can update the names of their own molecules"
        ):
            response = self.client.patch(
                reverse("molecule-detail", args=[self._normal_molecule.id]),
                {"name": "test"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_200_OK)
            self._normal_molecule.refresh_from_db()
            self.assertEqual(self._normal_molecule.name, "test")

        with self.subTest("Authenticated users cannot update other users' molecules"):
            response = self.client.put(
                reverse("molecule-detail", args=[self._admin_molecule.id]),
                {"name": "test", "smiles": "CC"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
            self.assertEqual(self._admin_molecule.name, "ethane")

    def test_molecule_admin_permissions(self) -> None:
        """Test that admin users can read and update all molecules."""
        self.client.force_authenticate(user=self._admin_user)

        with self.subTest("Admins are allowed to create molecules"):
            response = self.client.post(
                reverse("molecule-list"),
                {"name": "ethane", "smiles": "CC"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_201_CREATED)

        with self.subTest("Admins can update other users' molecules"):
            response = self.client.put(
                reverse("molecule-detail", args=[self._normal_molecule.id]),
                {"name": "test2", "smiles": "C"},
                format="json",
            )
            self.assertEqual(response.status_code, status.HTTP_200_OK)
            self._normal_molecule.refresh_from_db()
            self.assertEqual(self._normal_molecule.name, "test2")

    # def test_reaction_permissions(self) -> None:
    #     pass

    # def test_reaction_component_permissions(self) -> None:
    #     pass

    def test_user_profile_anonymous_permissions(self) -> None:
        """Test that anonymous users cannot read user profiles."""
        self.client.force_authenticate(user=self._anonymous_user)

        with self.subTest("Anonymous users cannot list user profiles"):
            response = self.client.get(reverse("userprofile-list"))
            self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

    def test_user_profile_authenticated_permissions(self) -> None:
        """Test that authenticated users can only read their own user profiles."""
        self.client.force_authenticate(user=self._normal_user)

        with self.subTest("Authenticated users can retrieve their own profile"):
            response = self.client.get(
                reverse("userprofile-detail", args=[self._normal_user.id])
            )
            self.assertEqual(response.status_code, status.HTTP_200_OK)

        with self.subTest("Authenticated users cannot retrieve other users' profiles"):
            response = self.client.get(
                reverse("userprofile-detail", args=[self._admin_user.id])
            )
            self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

        with self.subTest("Authenticated users cannot list profiles"):
            response = self.client.get(reverse("userprofile-list"))
            self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

    def test_user_profile_admin_permissions(self) -> None:
        """Test that admin users can read all user profiles."""
        self.client.force_authenticate(user=self._admin_user)

        with self.subTest("Admins can retrieve other users' profiles"):
            response = self.client.get(
                reverse("userprofile-detail", args=[self._normal_user.id])
            )
            self.assertEqual(response.status_code, status.HTTP_200_OK)

        with self.subTest("Admins can list profiles"):
            response = self.client.get(reverse("userprofile-list"))
            self.assertEqual(response.status_code, status.HTTP_200_OK)


# class UserAccountTest(APITestCase):
#     def setUp(self) -> None:
#         """Set up the dummy user for these tests."""
#         self._user = User.objects.create_user("test")

# def test_create_user(self):
#     pass
