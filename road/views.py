"""
The views and ViewSets for the ROAD REST API.
"""

# pylint: disable=too-many-ancestors

from __future__ import annotations

import logging
import os
from typing import NoReturn

from django.conf import settings
from django.contrib.auth.models import User  # pylint: disable=imported-auth-user
from django.db.models import Q, QuerySet
from django.db.utils import DataError
from query_parser import (  # pylint: disable=import-error, no-name-in-module
    QueryParserError,
    build_molecule_query,
)
from rest_framework.exceptions import NotFound
from rest_framework.generics import ListAPIView
from rest_framework.permissions import AllowAny
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.serializers import BaseSerializer
from rest_framework.views import APIView
from rest_framework.viewsets import ModelViewSet, ReadOnlyModelViewSet

from .access_policies import (
    MoleculeAccessPolicy,
    OverrideAccessViewSetMixin,
    ReactionAccessPolicy,
    ReactionComponentAccessPolicy,
    UserProfileAccessPolicy,
)
from .exceptions import InvalidQuery, ParameterNotProvided
from .models import Molecule, Reaction, ReactionComponent, UserProfile
from .permissions import ReadOnly
from .serializers import (
    MoleculeSerializer,
    ReactionComponentSerializer,
    ReactionSerializer,
    UserProfileSerializer,
)

logger = logging.getLogger(__name__)


class HideUnauthorisedMixin:  # pylint: disable=too-few-public-methods
    """Override the default permission_denied method to return a 404 instead of a 403."""

    def permission_denied(
        self, request: Request, message: str | None = None, code: str | None = None
    ) -> NoReturn:
        """
        Instead of returning a 403, return a 404, to conceal the existence of the resource.
        """
        raise NotFound()


class MoleculeViewSet(OverrideAccessViewSetMixin, ModelViewSet):  # type: ignore  # pylint: disable=too-few-public-methods
    """
    ViewSet for the Molecule model.
    """

    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer
    access_policy = MoleculeAccessPolicy

    def perform_create(self, serializer: BaseSerializer[Molecule]) -> None:
        """Create a new Molecule, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class ReactionViewSet(OverrideAccessViewSetMixin, ModelViewSet):  # type: ignore  # pylint: disable=too-few-public-methods
    """
    ViewSet for the Reaction model.
    """

    queryset = Reaction.objects.all()
    serializer_class = ReactionSerializer
    access_policy = ReactionAccessPolicy

    def perform_create(self, serializer: BaseSerializer[Reaction]) -> None:
        """Create a new Reaction, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class ReactionComponentViewSet(OverrideAccessViewSetMixin, ModelViewSet):  # type: ignore  # pylint: disable=too-few-public-methods
    """
    ViewSet for the ReactionComponent model.
    """

    queryset = ReactionComponent.objects.all()
    serializer_class = ReactionComponentSerializer
    access_policy = ReactionComponentAccessPolicy

    def perform_create(self, serializer: BaseSerializer[ReactionComponent]) -> None:
        """Create a new ReactionComponent, and set the owner to the request's user."""
        serializer.save(owner=self.request.user)


class UserProfileViewSet(
    OverrideAccessViewSetMixin,
    HideUnauthorisedMixin,
    ReadOnlyModelViewSet,  # type: ignore
):
    """
    ViewSet for the UserProfile model.
    """

    queryset = UserProfile.objects.all()
    serializer_class = UserProfileSerializer
    access_policy = UserProfileAccessPolicy


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

        queryset = Molecule.objects.filter(parsed_query)

        try:
            list(queryset)
        except DataError as error:
            raise InvalidQuery from error

        return queryset


class FlushView(APIView):
    """
    A view to flush the database for testing purposes.
    """

    permission_classes = [AllowAny]

    def post(self, request: Request) -> Response:
        """Flush the database if the user is allowed to."""
        if self.can_flush(request):
            self.flush_database()
            return Response(status=204)
        return Response(status=400)

    def can_flush(self, request: Request) -> bool:
        """Return True if the database can be flushed."""
        try:
            settings_allow_flush = settings.ALLOW_REMOTE_DATABASE_FLUSH
        except AttributeError:
            return False

        environment_allow_flush = (
            os.environ.get("ALLOW_REMOTE_DATABASE_FLUSH", "False") == "True"
        )

        try:
            secret = settings.REMOTE_DATABASE_FLUSH_SECRET
        except AttributeError:
            return False

        return (
            settings_allow_flush
            and environment_allow_flush
            and request.data.get("secret") == secret
        )

    def flush_database(self) -> None:
        """Flush the default database."""
        logger.warning("Flushing database")
        Reaction.objects.all().delete()
        Molecule.objects.all().delete()
        User.objects.all().delete()
        logger.warning("Database flushed")
