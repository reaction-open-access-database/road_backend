from django.db.models import Q
from rest_framework import viewsets
from rest_framework.exceptions import NotFound
from rest_framework.generics import ListAPIView
from query_parser import build_molecule_query, QueryParserError

from .exceptions import ParameterNotProvided, InvalidQuery
from .models import Molecule, Reaction, ReactionComponent, UserProfile
from .permissions import IsOwner, IsSuperUser, ReadOnly
from .serializers import (
    MoleculeSerializer,
    ReactionSerializer,
    ReactionComponentSerializer,
    UserProfileSerializer,
)


class HideUnauthorised:
    def permission_denied(self, request, message=None, code=None):
        raise NotFound()


class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    permission_classes = [IsSuperUser | IsOwner | ReadOnly]

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    serializer_class = ReactionSerializer
    permission_classes = [IsSuperUser | ReadOnly]

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)


class ReactionComponentViewSet(viewsets.ModelViewSet):
    queryset = ReactionComponent.objects.all()
    serializer_class = ReactionComponentSerializer
    permission_classes = [IsSuperUser | ReadOnly]

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)


class UserViewSet(HideUnauthorised, viewsets.ReadOnlyModelViewSet):
    queryset = UserProfile.objects.all()
    serializer_class = UserProfileSerializer

    def get_permissions(self):
        # Allow users to view and edit their own profiles,
        # but only superusers to view and edit all profiles
        if self.action == "retrieve":
            self.permission_classes = [IsOwner | IsSuperUser]
        elif self.action == "list":
            self.permission_classes = [IsSuperUser]
        return super().get_permissions()


class QueryView(ListAPIView):
    def _get_query(self):
        query = self.request.query_params.get("query")
        if not query:
            raise ParameterNotProvided("query parameter is required")

        return query


class MoleculeQueryView(QueryView):
    serializer_class = MoleculeSerializer
    permission_classes = [ReadOnly]

    def get_queryset(self):
        query = self._get_query()

        try:
            parsed_query = build_molecule_query(query, Q)
        except QueryParserError as error:
            raise InvalidQuery(error)

        return Molecule.objects.filter(parsed_query)
