"""helper_pta package

Public API:
  - Preprocess
  - ThermodynamicSpaceMod
  - ThermodynamicSpaceBasismod
  - SupTFSmodel
  - main

This package contains modularized components extracted from the original large
`helper_pta.py` script: preprocessing, thermodynamic-space helpers, basis
computation, the sampling TFS model and the orchestration `main` function.

Each module has a short module-level docstring and well-documented classes
and functions suitable for unit testing.
"""

from .preprocess import Preprocess
from .thermospace_mod import ThermodynamicSpaceMod
from .tbasis_mod import ThermodynamicSpaceBasismod
from .tfsmodel_mod import SupTFSmodel
from .main import main

__all__ = [
    "Preprocess",
    "ThermodynamicSpaceMod",
    "ThermodynamicSpaceBasismod",
    "SupTFSmodel",
    "main",
]
