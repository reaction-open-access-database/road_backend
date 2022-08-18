from rest_framework import viewsets, permissions
from .models import Molecule
from .serializers import MoleculeSerializer


class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    permission_classes = [permissions.IsAuthenticated]
