"""
Trinarization methods for converting continuous data to ternary values.
"""

import numpy as np
from sklearn.cluster import KMeans
from .results import TrinarizationResult


def trinarize_kmeans(
    vect,
    nstart=1,
    iter_max=10,
    na_rm=False,
    random_state=None
):
    """
    Trinarize a real-valued vector using k-Means clustering.
    
    This method uses k-Means clustering with k=3 to separate the data into
    three clusters, then assigns ternary values (0, 1, and 2) based on cluster centers.
    
    Parameters
    ----------
    vect : array-like
        Input vector of numerical values to trinarize
    nstart : int, optional (default=1)
        Number of random initializations for k-Means (n_init parameter)
    iter_max : int, optional (default=10)
        Maximum number of iterations for k-Means
    na_rm : bool, optional (default=False)
        If True, remove NA/NaN values before trinarization
    random_state : int, optional (default=None)
        Random seed for reproducibility
        
    Returns
    -------
    TrinarizationResult
        Object containing:
        - original_measurements: original data values
        - trinarized_measurements: ternary values (0, 1, or 2)
        - threshold1: first threshold (separating 0 from 1)
        - threshold2: second threshold (separating 1 from 2)
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
    >>> from bitrina import trinarize_kmeans
    >>> data = np.array([1.1, 1.2, 3.5, 3.6, 6.1, 6.2])
    >>> result = trinarize_kmeans(data)
    >>> print(result.trinarized_measurements)
    [0 0 1 1 2 2]
    >>> print(result.threshold1, result.threshold2)
    2.35 4.8
    """
    # Convert to numpy array
    vect = np.asarray(vect, dtype=float)
    
    # Check for NA values
    if na_rm:
        vect = vect[~np.isnan(vect)]
    elif np.any(np.isnan(vect)):
        raise ValueError("Cannot trinarize in the presence of NA values!")
    
    # Check for Inf values
    if np.any(~np.isfinite(vect)):
        raise ValueError("Cannot trinarize Inf values!")
    
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
    
    # Perform k-Means clustering with k=3
    kmeans = KMeans(
        n_clusters=3,
        n_init=nstart,
        max_iter=iter_max,
        random_state=random_state
    )
    
    # Reshape for sklearn (requires 2D input)
    vect_2d = vect.reshape(-1, 1)
    cluster_labels = kmeans.fit_predict(vect_2d)
    centers = kmeans.cluster_centers_.flatten()
    
    # Sort centers and reassign labels accordingly
    center_order = np.argsort(centers)
    trinarized = np.zeros_like(cluster_labels)
    
    for new_label, old_label in enumerate(center_order):
        trinarized[cluster_labels == old_label] = new_label
    
    # Calculate thresholds
    low_values = vect[trinarized == 0]
    mid_values = vect[trinarized == 1]
    high_values = vect[trinarized == 2]
    
    threshold1 = (np.max(low_values) + np.min(mid_values)) / 2
    threshold2 = (np.max(mid_values) + np.min(high_values)) / 2
    
    # Return result object
    return TrinarizationResult(
        original_measurements=vect,
        trinarized_measurements=trinarized.astype(int),
        threshold1=threshold1,
        threshold2=threshold2,
        p_value=None,  # Dip test not implemented in basic version
        method="k-Means"
    )
