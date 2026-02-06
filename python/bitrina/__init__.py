"""
BiTrinA - Binarization and Trinarization of One-Dimensional Data

A Python implementation of the BiTrinA package for binarization and trinarization
of continuous data, with special support for single-cell RNA-seq data (AnnData objects).
"""

from .binarize import binarize_kmeans
from .trinarize import trinarize_kmeans
from .matrix import binarize_matrix, trinarize_matrix
from .results import BinarizationResult, TrinarizationResult

__version__ = "0.1.0"
__all__ = [
    "binarize_kmeans",
    "trinarize_kmeans",
    "binarize_matrix",
    "trinarize_matrix",
    "BinarizationResult",
    "TrinarizationResult",
]
