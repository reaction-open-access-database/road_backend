"""
Wraps some RDKit functions to raise an exception if there is an error.
"""

import logging
import re
import sys
from contextlib import contextmanager
from io import StringIO
from typing import Iterator, List, Optional

from rdkit.Chem import Mol, MolFromSmiles
from rdkit.rdBase import LogToPythonStderr

from .exceptions import InvalidMolecule

LogToPythonStderr()
logger = logging.getLogger(__name__)

__all__ = [
    "smiles_to_mol",
]


def smiles_to_mol(smiles: str) -> Mol:
    """
    Converts a SMILES string to an RDKit molecule.
    Raises an exception if the conversion fails.
    """
    errors: List[str] = []

    with get_rdkit_error_lines(errors):
        mol = MolFromSmiles(smiles)

    if mol is None:
        if len(errors) == 0:
            logger.warning("Invalid SMILES: {} (no error messages returned)", smiles)
            raise InvalidMolecule("Invalid SMILES")
        raise InvalidMolecule(errors)

    if len(errors) > 0:
        logger.warning("Valid SMILES: {} (with error messages: {})", smiles, errors)
        raise InvalidMolecule(errors)

    return mol


@contextmanager
def get_rdkit_error_lines(errors: List[str]) -> Iterator[None]:
    """
    Context manager that captures error messages, and splits them into lines,
    which are appended to the errors list.
    """

    original = sys.stderr

    error_io = sys.stderr = StringIO()

    try:
        yield
    finally:
        sys.stderr = original
        for line in error_io.getvalue().splitlines():
            error = parse_rdkit_error_line(line)
            if error is not None:
                errors.append(error)


def parse_rdkit_error_line(line: str) -> Optional[str]:
    """
    Parses an RDKit error line, and returns the error message.
    Returns None if the line has an unexpected format.

    Example: "[16:04:45] Explicit valence for atom # 7 N, 4, is greater than permitted"
    becomes "Explicit valence for atom # 7 N, 4, is greater than permitted"
    """
    match = re.fullmatch(r"\[\d\d:\d\d:\d\d] (.*)", line)
    if match is None:
        return None
    return match.group(1)
