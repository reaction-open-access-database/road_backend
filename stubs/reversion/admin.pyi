from typing import Any

from django.contrib.admin import ModelAdmin

class VersionAdmin(ModelAdmin):  # type: ignore
    pass

def __getattr__(name: str) -> Any: ...
