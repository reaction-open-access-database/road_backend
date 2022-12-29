"""
Defines the admin interface for ROAD.
"""

from django.contrib import admin
from .models import UserProfile, Molecule, Reaction, ReactionComponent, ReactionSource


admin.site.register(UserProfile)

admin.site.register(Reaction)
admin.site.register(ReactionComponent)
admin.site.register(ReactionSource)


@admin.register(Molecule)
class MoleculeAdmin(admin.ModelAdmin[Molecule]):
    """
    Custom admin interface for the Molecule model.
    Displays the molecule name and molecular formula.
    """

    list_display = ("name", "molecular_formula")
