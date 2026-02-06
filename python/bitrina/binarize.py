"""
Binarization methods for converting continuous data to binary values.
"""

import numpy as np
from sklearn.cluster import KMeans
from .results import BinarizationResult


def binarize_kmeans(
    vect,
    nstart=1,
    iter_max=10,
    na_rm=False,
    random_state=None
):
    """
    Binarize a real-valued vector using k-Means clustering.
    
    This method uses k-Means clustering with k=2 to separate the data into
    two clusters, then assigns binary values (0 and 1) based on cluster centers.
    
    Parameters
    ----------
    vect : array-like
        Input vector of numerical values to binarize
    nstart : int, optional (default=1)
        Number of random initializations for k-Means (n_init parameter)
    iter_max : int, optional (default=10)
        Maximum number of iterations for k-Means
    na_rm : bool, optional (default=False)
        If True, remove NA/NaN values before binarization
    random_state : int, optional (default=None)
        Random seed for reproducibility
        
    Returns
    -------
    BinarizationResult
        Object containing:
        - original_measurements: original data values
        - binarized_measurements: binary values (0 or 1)
        - threshold: computed threshold between clusters
        - p_value: None (statistical test not implemented)
        - method: "k-Means"
        
    Raises
    ------
    ValueError
        If input is not numeric, constant, has too few elements,
        contains NA values (when na_rm=False), or contains Inf values
        
    Examples
    --------
    >>> import numpy as np
    >>> from bitrina import binarize_kmeans
    >>> data = np.array([1.2, 1.5, 1.3, 5.6, 5.9, 5.4])
    >>> result = binarize_kmeans(data)
    >>> print(result.binarized_measurements)
    [0 0 0 1 1 1]
    >>> print(result.threshold)
    3.45
    """
    # Convert to numpy array
    vect = np.asarray(vect, dtype=float)
    
    # Check for NA values
    if na_rm:
        vect = vect[~np.isnan(vect)]
    elif np.any(np.isnan(vect)):
        raise ValueError("Cannot binarize in the presence of NA values!")
    
    # Check for Inf values
    if np.any(~np.isfinite(vect)):
        raise ValueError("Cannot binarize Inf values!")
    
    # Check vector length
    if len(vect) < 3:
        raise ValueError("The input vector must have at least 3 entries!")
    
    # Check if constant
    if len(np.unique(vect)) == 1:
        raise ValueError("The input vector is constant!")
    
    # Validate parameters
    if not isinstance(nstart, int) or nstart < 0:
        raise ValueError("'nstart' must be a non-negative integer!")
    if not isinstance(iter_max, int) or iter_max < 0:
        raise ValueError("'iter_max' must be a non-negative integer!")
    
    # Perform k-Means clustering with k=2
    kmeans = KMeans(
        n_clusters=2,
        n_init=nstart,
        max_iter=iter_max,
        random_state=random_state
    )
    
    # Reshape for sklearn (requires 2D input)
    vect_2d = vect.reshape(-1, 1)
    cluster_labels = kmeans.fit_predict(vect_2d)
    centers = kmeans.cluster_centers_.flatten()
    
    # Assign binary values: higher center gets 1, lower gets 0
    if centers[0] > centers[1]:
        binarized = np.abs(cluster_labels - 1)  # Flip labels
    else:
        binarized = cluster_labels
    
    # Calculate threshold as midpoint between min of high cluster and max of low cluster
    low_values = vect[binarized == 0]
    high_values = vect[binarized == 1]
    threshold = (np.max(low_values) + np.min(high_values)) / 2
    
    # Return result object
    return BinarizationResult(
        original_measurements=vect,
        binarized_measurements=binarized.astype(int),
        threshold=threshold,
        p_value=None,  # Dip test not implemented in basic version
        method="k-Means"
    )
