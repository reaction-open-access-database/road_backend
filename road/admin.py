from django.contrib import admin
from .models import UserProfile, Molecule, Reaction, ReactionComponent, \
    ReactionSource
from reversion.admin import VersionAdmin


admin.site.register(UserProfile)


@admin.register(Reaction)
class ReactionAdmin(VersionAdmin):
    pass


@admin.register(ReactionComponent)
class ReactionComponentAdmin(VersionAdmin):
    pass


@admin.register(ReactionSource)
class ReactionSourceAdmin(VersionAdmin):
    pass


@admin.register(Molecule)
class MoleculeAdmin(VersionAdmin):
    list_display = ('name', 'molecular_formula')
