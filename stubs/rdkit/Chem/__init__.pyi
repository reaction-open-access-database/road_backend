from typing import Any

from rdkit.AllChem import Mol

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
def __getattr__(name: str) -> Any: ...
