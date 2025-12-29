"""ChEMBL-PDB Linker: Link bioactivity data with structural information."""

__version__ = "0.1.0"

from chembl_pdb_linker.config import Config
from chembl_pdb_linker.pipeline import Pipeline

__all__ = ["Config", "Pipeline", "__version__"]
