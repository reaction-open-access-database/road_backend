from typing import Any

from rdkit.Chem.rdchem import Mol

from . import rdMolDescriptors

def MolToInchi(
    mol: Mol, options: str = ..., logLevel: Any = ..., treatWarningAsError: bool = ...
) -> str: ...
def MolToSmiles(
    mol: Mol,
    isomericSmiles: bool = ...,
    kekuleSmiles: bool = ...,
    rootedAtAtom: int = ...,
    canonical: bool = ...,
    allBondsExplicit: bool = ...,
    allHsExplicit: bool = ...,
    doRandom: bool = ...,
) -> str: ...
def MolToJSON(mol: Mol) -> str: ...
def __getattr__(name: str) -> Any: ...

__all__ = [
    "MolToInchi",
    "MolToSmiles",
    "MolToJSON",
    "Mol",
]
