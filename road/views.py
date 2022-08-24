from rest_framework import viewsets
from .models import Molecule
from .serializers import MoleculeSerializer
from .permissions import IsOwnerOrReadOnly


class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    permission_classes = [IsOwnerOrReadOnly]
