from typing import Any, Optional

from rdkit.Chem.AllChem import AtomPairsParameters, Mol

class MolDraw2D:
    def DrawMolecule(
        self,
        mol: Mol,
        highlightAtoms: Optional[AtomPairsParameters] = ...,
        highlightAtomColors: Optional[AtomPairsParameters] = ...,
        highlightAtomRadii: Optional[AtomPairsParameters] = ...,
        confId: int = ...,
        legend: str = ...,
    ) -> None: ...

class MolDraw2DSVG(MolDraw2D):
    def __init__(
        self,
        width: int,
        height: int,
        panelWidth: int = ...,
        panelHeight: int = ...,
        noFreetype: bool = ...,
    ) -> None: ...
    def FinishDrawing(self) -> None: ...
    def GetDrawingText(self) -> str: ...

def __getattr__(name: str) -> Any: ...
