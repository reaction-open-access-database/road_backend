"""
Defines the admin interface for ROAD.
"""

import json
from typing import Any

from django.contrib import admin
from django.forms import JSONField, ModelForm
from reversion.admin import VersionAdmin

from .models import Molecule, Reaction, ReactionComponent, ReactionSource, UserProfile

admin.site.register(UserProfile)


@admin.register(Reaction)
class ReactionAdmin(VersionAdmin):
    """Admin interface for the Reaction model."""


@admin.register(ReactionComponent)
class ReactionComponentAdmin(VersionAdmin):
    """Admin interface for the ReactionComponent model."""


@admin.register(ReactionSource)
class ReactionSourceAdmin(VersionAdmin):
    """Admin interface for the ReactionSource model."""


class PrettyJSONEncoder(json.JSONEncoder):
    """JSON encoder that formats the JSON nicely (indented and sorted)."""

    def __init__(self, indent: Any, sort_keys: Any, **kwargs: Any) -> None:
        super().__init__(indent=2, sort_keys=True, **kwargs)


class MoleculeForm(
    ModelForm[Molecule]
):  # pylint: disable=unsubscriptable-object, too-few-public-methods
    """Form for the Molecule model."""

    molecule = JSONField(encoder=PrettyJSONEncoder)


@admin.register(Molecule)
class MoleculeAdmin(
    VersionAdmin
):  # pylint: disable=unsubscriptable-object, too-few-public-methods
    """
    Custom admin interface for the Molecule model.
    Displays the molecule name and molecular formula.
    """

    list_display = ("id", "name", "molecular_formula")
    form = MoleculeForm
    readonly_fields = ("molecular_formula",)
