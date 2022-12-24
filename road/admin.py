from django.contrib import admin
from .models import UserProfile, Molecule, Reaction, ReactionComponent, ReactionSource


admin.site.register(UserProfile)

admin.site.register(Reaction)
admin.site.register(ReactionComponent)
admin.site.register(ReactionSource)


@admin.register(Molecule)
class MoleculeAdmin(admin.ModelAdmin):
    list_display = ("name", "molecular_formula")
