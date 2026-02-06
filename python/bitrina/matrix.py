"""
Matrix operations for batch binarization and trinarization.

Supports processing multiple vectors at once, with special support for
single-cell RNA-seq data stored in AnnData objects.
"""

import numpy as np
import pandas as pd
from typing import Union, Optional, Literal
from .binarize import binarize_kmeans
from .trinarize import trinarize_kmeans


def binarize_matrix(
    mat: Union[np.ndarray, pd.DataFrame, 'anndata.AnnData'],
    method: Literal["kMeans"] = "kMeans",
    layer: Optional[str] = None,
    **kwargs
):
    """
    Binarize a matrix where each row is processed independently.
    
    This function applies binarization to each row of a matrix independently,
    which is useful for processing gene expression data where rows represent
    genes/features and columns represent samples/cells.
    
    Parameters
    ----------
    mat : np.ndarray, pd.DataFrame, or AnnData
        Input matrix to binarize. If AnnData, will use .X by default or
        specified layer. Rows are features, columns are samples.
    method : {"kMeans"}, optional (default="kMeans")
        Binarization method to use. Currently only "kMeans" is implemented.
    layer : str, optional (default=None)
        For AnnData input, which layer to use. If None, uses .X
    **kwargs
        Additional arguments passed to the binarization function
        
    Returns
    -------
    pd.DataFrame
        DataFrame with binarized values for each row, plus columns for:
        - threshold: computed threshold for each row
        - p_value: p-value for each row (if applicable)
        
    Examples
    --------
    >>> import numpy as np
    >>> from bitrina import binarize_matrix
    >>> # Matrix with 3 genes (rows) and 6 cells (columns)
    >>> mat = np.array([
    ...     [1.2, 1.5, 1.3, 5.6, 5.9, 5.4],
    ...     [0.1, 0.2, 0.15, 2.1, 2.3, 2.0],
    ...     [3.1, 3.2, 3.0, 7.5, 7.8, 7.6]
    ... ])
    >>> result = binarize_matrix(mat)
    >>> print(result.iloc[:, :6])  # First 6 columns (binarized values)
       V0  V1  V2  V3  V4  V5
    0   0   0   0   1   1   1
    1   0   0   0   1   1   1
    2   0   0   0   1   1   1
    """
    # Handle AnnData input
    if hasattr(mat, 'X'):  # AnnData-like object
        import anndata
        if not isinstance(mat, anndata.AnnData):
            raise ValueError("Input appears to be AnnData-like but is not an AnnData object")
        
        if layer is not None:
            if layer not in mat.layers:
                raise ValueError(f"Layer '{layer}' not found in AnnData object")
            data = mat.layers[layer]
        else:
            data = mat.X
            
        # Convert sparse to dense if needed
        if hasattr(data, 'toarray'):
            data = data.toarray()
        
        # Store var_names (gene names) for output
        row_names = mat.var_names if mat.var_names is not None else None
        col_names = mat.obs_names if mat.obs_names is not None else None
    else:
        # Handle numpy array or pandas DataFrame
        if isinstance(mat, pd.DataFrame):
            data = mat.values
            row_names = mat.index
            col_names = mat.columns
        else:
            data = np.asarray(mat)
            row_names = None
            col_names = None
    
    # Select binarization function
    if method == "kMeans":
        bin_func = binarize_kmeans
    else:
        raise ValueError(f"Unknown method: {method}. Currently only 'kMeans' is supported.")
    
    # Process each row
    results = []
    for i, row in enumerate(data):
        try:
            result = bin_func(row, **kwargs)
            row_result = list(result.binarized_measurements) + [
                result.threshold,
                result.p_value
            ]
            results.append(row_result)
        except ValueError as e:
            # If binarization fails for a row, fill with NaN
            n_cols = data.shape[1]
            row_result = [np.nan] * n_cols + [np.nan, np.nan]
            results.append(row_result)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Set column names
    if col_names is not None:
        df.columns = list(col_names) + ['threshold', 'p_value']
    else:
        n_cols = data.shape[1]
        df.columns = [f'V{i}' for i in range(n_cols)] + ['threshold', 'p_value']
    
    # Set row names
    if row_names is not None:
        df.index = row_names
    
    return df


def trinarize_matrix(
    mat: Union[np.ndarray, pd.DataFrame, 'anndata.AnnData'],
    method: Literal["kMeans"] = "kMeans",
    layer: Optional[str] = None,
    **kwargs
):
    """
    Trinarize a matrix where each row is processed independently.
    
    This function applies trinarization to each row of a matrix independently,
    which is useful for processing gene expression data where rows represent
    genes/features and columns represent samples/cells.
    
    Parameters
    ----------
    mat : np.ndarray, pd.DataFrame, or AnnData
        Input matrix to trinarize. If AnnData, will use .X by default or
        specified layer. Rows are features, columns are samples.
    method : {"kMeans"}, optional (default="kMeans")
        Trinarization method to use. Currently only "kMeans" is implemented.
    layer : str, optional (default=None)
        For AnnData input, which layer to use. If None, uses .X
    **kwargs
        Additional arguments passed to the trinarization function
        
    Returns
    -------
    pd.DataFrame
        DataFrame with trinarized values for each row, plus columns for:
        - threshold1: first threshold for each row
        - threshold2: second threshold for each row
        - p_value: p-value for each row (if applicable)
        
    Examples
    --------
    >>> import numpy as np
    >>> from bitrina import trinarize_matrix
    >>> # Matrix with 2 genes (rows) and 6 cells (columns)
    >>> mat = np.array([
    ...     [1.1, 1.2, 3.5, 3.6, 6.1, 6.2],
    ...     [0.5, 0.6, 2.5, 2.6, 5.5, 5.6]
    ... ])
    >>> result = trinarize_matrix(mat)
    >>> print(result.iloc[:, :6])  # First 6 columns (trinarized values)
       V0  V1  V2  V3  V4  V5
    0   0   0   1   1   2   2
    1   0   0   1   1   2   2
    """
    # Handle AnnData input
    if hasattr(mat, 'X'):  # AnnData-like object
        import anndata
        if not isinstance(mat, anndata.AnnData):
            raise ValueError("Input appears to be AnnData-like but is not an AnnData object")
        
        if layer is not None:
            if layer not in mat.layers:
                raise ValueError(f"Layer '{layer}' not found in AnnData object")
            data = mat.layers[layer]
        else:
            data = mat.X
            
        # Convert sparse to dense if needed
        if hasattr(data, 'toarray'):
            data = data.toarray()
        
        # Store var_names (gene names) for output
        row_names = mat.var_names if mat.var_names is not None else None
        col_names = mat.obs_names if mat.obs_names is not None else None
    else:
        # Handle numpy array or pandas DataFrame
        if isinstance(mat, pd.DataFrame):
            data = mat.values
            row_names = mat.index
            col_names = mat.columns
        else:
            data = np.asarray(mat)
            row_names = None
            col_names = None
    
    # Select trinarization function
    if method == "kMeans":
        tri_func = trinarize_kmeans
    else:
        raise ValueError(f"Unknown method: {method}. Currently only 'kMeans' is supported.")
    
    # Process each row
    results = []
    for i, row in enumerate(data):
        try:
            result = tri_func(row, **kwargs)
            row_result = list(result.trinarized_measurements) + [
                result.threshold1,
                result.threshold2,
                result.p_value
            ]
            results.append(row_result)
        except ValueError as e:
            # If trinarization fails for a row, fill with NaN
            n_cols = data.shape[1]
            row_result = [np.nan] * n_cols + [np.nan, np.nan, np.nan]
            results.append(row_result)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Set column names
    if col_names is not None:
        df.columns = list(col_names) + ['threshold1', 'threshold2', 'p_value']
    else:
        n_cols = data.shape[1]
        df.columns = [f'V{i}' for i in range(n_cols)] + ['threshold1', 'threshold2', 'p_value']
    
    # Set row names
    if row_names is not None:
        df.index = row_names
    
    return df
