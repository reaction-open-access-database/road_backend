"""
Defines the admin interface for ROAD.
"""

from django.contrib import admin
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


@admin.register(Molecule)
class MoleculeAdmin(
    VersionAdmin
):  # pylint: disable=unsubscriptable-object, too-few-public-methods
    """
    Custom admin interface for the Molecule model.
    Displays the molecule name and molecular formula.
    """

    list_display = ("name", "molecular_formula")
