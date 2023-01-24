"""
Wraps some RDKit functions to raise an exception if there is an error.
"""

import logging
import re
import sys
from contextlib import contextmanager
from io import StringIO
from typing import Callable, Iterator, List, Optional

from rdkit.Chem import Mol, MolFromInchi, MolFromSmiles
from rdkit.rdBase import LogToPythonStderr

from .exceptions import InvalidMolecule

LogToPythonStderr()
logger = logging.getLogger(__name__)

__all__ = [
    "smiles_to_mol",
    "inchi_to_mol",
]


def smiles_to_mol(smiles: str) -> Mol:
    """
    Converts a SMILES string to an RDKit molecule.
    Raises an exception if the conversion fails.
    """
    return to_mol(MolFromSmiles, "SMILES", smiles)


def inchi_to_mol(inchi: str) -> Mol:
    """
    Converts an InChI string to an RDKit molecule.
    Raises an exception if the conversion fails.
    """
    return to_mol(MolFromInchi, "InChI", inchi)


def to_mol(
    func: Callable[[str], Optional[Mol]], input_type: str, input_value: str
) -> Mol:
    """
    Convert into an RDKit molecule.
    Raises an exception if the conversion fails.
    """
    errors: List[str] = []

    with get_rdkit_error_lines(errors):
        mol = func(input_value)

    # When parsing invalid InChI strings, RDKit returns the "ERROR: " string, which isn't helpful
    # to the user. We'll remove it from the error messages.
    errors.remove("ERROR: ")

    if mol is None:
        if len(errors) == 0:
            logger.warning(
                "Invalid %s: %s (no error messages returned)", input_type, input_value
            )
            raise InvalidMolecule(f"Invalid {input_type}")
        raise InvalidMolecule(errors)

    if len(errors) > 0:
        logger.warning(
            "Valid %s: %s (with error messages: %s)", input_type, input_value, errors
        )
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
