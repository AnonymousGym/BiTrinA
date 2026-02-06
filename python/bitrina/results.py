"""
Result classes for storing binarization and trinarization outputs.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class BinarizationResult:
    """
    Result object for binarization operations.
    
    Attributes
    ----------
    original_measurements : np.ndarray
        The original continuous data values
    binarized_measurements : np.ndarray
        The binarized values (0 or 1)
    threshold : float
        The threshold used for binarization
    p_value : Optional[float]
        P-value from statistical test (if applicable)
    method : str
        The method used for binarization
    """
    original_measurements: np.ndarray
    binarized_measurements: np.ndarray
    threshold: float
    p_value: Optional[float] = None
    method: str = "unknown"
    
    def __repr__(self):
        return (f"BinarizationResult(method='{self.method}', "
                f"threshold={self.threshold:.4f}, "
                f"n_samples={len(self.original_measurements)})")


@dataclass
class TrinarizationResult:
    """
    Result object for trinarization operations.
    
    Attributes
    ----------
    original_measurements : np.ndarray
        The original continuous data values
    trinarized_measurements : np.ndarray
        The trinarized values (0, 1, or 2)
    threshold1 : float
        The first threshold (separating 0 from 1)
    threshold2 : float
        The second threshold (separating 1 from 2)
    p_value : Optional[float]
        P-value from statistical test (if applicable)
    method : str
        The method used for trinarization
    """
    original_measurements: np.ndarray
    trinarized_measurements: np.ndarray
    threshold1: float
    threshold2: float
    p_value: Optional[float] = None
    method: str = "unknown"
    
    def __repr__(self):
        return (f"TrinarizationResult(method='{self.method}', "
                f"threshold1={self.threshold1:.4f}, "
                f"threshold2={self.threshold2:.4f}, "
                f"n_samples={len(self.original_measurements)})")
