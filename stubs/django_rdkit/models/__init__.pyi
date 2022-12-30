from typing import Any

from django.db.models import *

class MolField:
    pass

def __getattr__(name: str) -> Any: ...
