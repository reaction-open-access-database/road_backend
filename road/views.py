"""
The views and ViewSets for the ROAD REST API.
"""

# pylint: disable=too-many-ancestors

from __future__ import annotations
from typing import List, NoReturn, Any

from django.db.models import Q, QuerySet
from rest_framework import viewsets
from rest_framework.exceptions import NotFound
from rest_framework.generics import ListAPIView
from rest_framework.request import Request
from rest_framework.serializers import BaseSerializer
from query_parser import (  # pylint: disable=import-error, no-name-in-module
    build_molecule_query,
    QueryParserError,
)

from .exceptions import ParameterNotProvided, InvalidQuery
from .models import Molecule, Reaction, ReactionComponent, UserProfile
from .permissions import IsOwner, IsSuperUser, ReadOnly
from .serializers import (
    MoleculeSerializer,
    ReactionSerializer,
    ReactionComponentSerializer,
    UserProfileSerializer,
)


class HideUnauthorised:  # pylint: disable=too-few-public-methods
    """Override the default permission_denied method to return a 404 instead of a 403."""

    def permission_denied(
        self, request: Request, message: str | None = None, code: str | None = None
    ) -> NoReturn:
        """
        Instead of returning a 403, return a 404, to conceal the existence of the resource.
        """
        raise NotFound()


class MoleculeViewSet(
    viewsets.ModelViewSet[Molecule]
):  # pylint: disable=too-few-public-methods
    """
    ViewSet for the Molecule model.
    """

    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    permission_classes = [IsSuperUser | IsOwner | ReadOnly]

    def perform_create(self, serializer: BaseSerializer[Molecule]) -> None:
        """Create a new Molecule, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class ReactionViewSet(
    viewsets.ModelViewSet[Reaction]
):  # pylint: disable=too-few-public-methods
    """
    ViewSet for the Reaction model.
    """

    queryset = Reaction.objects.all()
    serializer_class = ReactionSerializer
    permission_classes = [IsSuperUser | ReadOnly]

    def perform_create(self, serializer: BaseSerializer[Reaction]) -> None:
        """Create a new Reaction, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class ReactionComponentViewSet(
    viewsets.ModelViewSet[ReactionComponent]
):  # pylint: disable=too-few-public-methods
    """
    ViewSet for the ReactionComponent model.
    """

    queryset = ReactionComponent.objects.all()
    serializer_class = ReactionComponentSerializer
    permission_classes = [IsSuperUser | ReadOnly]

    def perform_create(self, serializer: BaseSerializer[ReactionComponent]) -> None:
        """Create a new ReactionComponent, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class UserViewSet(HideUnauthorised, viewsets.ReadOnlyModelViewSet[UserProfile]):
    """
    ViewSet for the UserProfile model.
    """

    queryset = UserProfile.objects.all()
    serializer_class = UserProfileSerializer
    permission_classes = [IsSuperUser]

    def get_permissions(self) -> List[Any]:
        """
        Allow users to view and edit their own profiles,
        but only superusers to view and edit all profiles.
        """
        if self.action == "retrieve":
            self.permission_classes = [IsOwner | IsSuperUser]
        elif self.action == "list":
            self.permission_classes = [IsSuperUser]
        return super().get_permissions()


class QueryView(ListAPIView):  # type: ignore
    """
    A view for querying the database.
    """

    def _get_query(self) -> str:
        """
        Ensure that the query parameter is provided.
        Returns the query parameter if present, otherwise raises ParameterNotProvided.
        """
        query = self.request.query_params.get("query")
        if not query:
            raise ParameterNotProvided("query parameter is required")

        return query


class MoleculeQueryView(QueryView):
    """
    A view for querying the Molecule model.
    """

    serializer_class = MoleculeSerializer
    permission_classes = [ReadOnly]

    def get_queryset(self) -> QuerySet[Molecule]:
        query = self._get_query()

        try:
            parsed_query = build_molecule_query(query, Q)
        except QueryParserError as error:
            raise InvalidQuery from error

        return Molecule.objects.filter(parsed_query)
