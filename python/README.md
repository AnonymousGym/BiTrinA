# BiTrinA - Python Implementation

A Python implementation of the BiTrinA (Binarization and Trinarization) package for converting continuous data into binary or ternary discrete values. This implementation is specifically designed to work seamlessly with single-cell RNA-seq data stored in AnnData objects.

## Overview

BiTrinA provides methods to discretize continuous measurements into:
- **Binarization**: Convert to two states (0 and 1)
- **Trinarization**: Convert to three states (0, 1, and 2)

This is particularly useful in bioinformatics for:
- Gene expression analysis
- Single-cell RNA-seq data preprocessing
- Boolean network inference
- Feature discretization for machine learning

## Installation

```bash
cd python
pip install -e .
```

For single-cell RNA-seq support with AnnData:
```bash
pip install -e ".[anndata]"
```

## Quick Start

### Basic Usage with NumPy Arrays

```python
import numpy as np
from bitrina import binarize_kmeans, trinarize_kmeans

# Binarization example
data = np.array([1.2, 1.5, 1.3, 5.6, 5.9, 5.4])
result = binarize_kmeans(data)
print(result.binarized_measurements)  # [0 0 0 1 1 1]
print(result.threshold)  # 3.45

# Trinarization example
data = np.array([1.1, 1.2, 3.5, 3.6, 6.1, 6.2])
result = trinarize_kmeans(data)
print(result.trinarized_measurements)  # [0 0 1 1 2 2]
print(result.threshold1, result.threshold2)  # 2.35 4.8
```

### Working with Matrices (Gene Expression Data)

```python
import numpy as np
from bitrina import binarize_matrix, trinarize_matrix

# Example: 3 genes (rows) × 6 cells (columns)
expression_data = np.array([
    [1.2, 1.5, 1.3, 5.6, 5.9, 5.4],  # Gene 1
    [0.1, 0.2, 0.15, 2.1, 2.3, 2.0],  # Gene 2
    [3.1, 3.2, 3.0, 7.5, 7.8, 7.6]   # Gene 3
])

# Binarize each gene independently
binarized = binarize_matrix(expression_data)
print(binarized)
```

### Single-Cell RNA-seq with AnnData

```python
import anndata as ad
import numpy as np
from bitrina import binarize_matrix, trinarize_matrix

# Load your single-cell data
adata = ad.read_h5ad("your_data.h5ad")

# Binarize gene expression across cells
# Each gene (row in adata.X transposed) is binarized independently
binarized = binarize_matrix(adata.T, method="kMeans")

# Or trinarize for three expression levels
trinarized = trinarize_matrix(adata.T, method="kMeans")

# Use a specific layer
binarized = binarize_matrix(adata.T, method="kMeans", layer="log1p")
```

## Input and Output Specifications

### Input Formats

BiTrinA accepts multiple input formats:

1. **Vectors** (for `binarize_kmeans`, `trinarize_kmeans`):
   - NumPy arrays: `np.array([1.2, 3.4, 5.6])`
   - Python lists: `[1.2, 3.4, 5.6]`
   - Pandas Series: `pd.Series([1.2, 3.4, 5.6])`

2. **Matrices** (for `binarize_matrix`, `trinarize_matrix`):
   - NumPy 2D arrays: `np.array([[1, 2], [3, 4]])`
   - Pandas DataFrames: `pd.DataFrame([[1, 2], [3, 4]])`
   - AnnData objects: `adata` (from scanpy/anndata)

### Output Formats

1. **For Vector Functions** (`binarize_kmeans`, `trinarize_kmeans`):
   - Returns `BinarizationResult` or `TrinarizationResult` objects with attributes:
     - `original_measurements`: Input data
     - `binarized_measurements` / `trinarized_measurements`: Discretized values
     - `threshold` (binarization) or `threshold1`, `threshold2` (trinarization)
     - `p_value`: Statistical test result (if applicable)
     - `method`: Method name used

2. **For Matrix Functions** (`binarize_matrix`, `trinarize_matrix`):
   - Returns `pandas.DataFrame` with:
     - Discretized values for each cell/sample (columns)
     - `threshold` or `threshold1`, `threshold2` columns
     - `p_value` column
     - Original row names (gene names) and column names (cell names) preserved

### Working with scRNA-seq AnnData

When working with AnnData objects from scanpy:

```python
import scanpy as sc
import anndata as ad
from bitrina import binarize_matrix

# Load standard scRNA-seq data
adata = sc.read_10x_h5("filtered_gene_bc_matrices_h5.h5")

# Preprocess
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Binarize gene expression (note: transpose to get genes as rows)
binarized_df = binarize_matrix(adata.T, method="kMeans")

# Result is a DataFrame with:
# - Rows: genes (from adata.var_names)
# - Columns: cells (from adata.obs_names) + threshold + p_value

# You can store back in AnnData
adata.layers["binarized"] = binarized_df.iloc[:, :-2].T.values
```

## Available Methods

Currently implemented:
- **k-Means clustering** (`method="kMeans"`): Fast and simple clustering-based approach

## API Reference

### Vector Functions

#### `binarize_kmeans(vect, nstart=1, iter_max=10, na_rm=False, random_state=None)`

Binarize a vector using k-Means clustering.

**Parameters:**
- `vect`: Input vector (array-like)
- `nstart`: Number of random initializations (default: 1)
- `iter_max`: Maximum iterations (default: 10)
- `na_rm`: Remove NA values (default: False)
- `random_state`: Random seed for reproducibility (default: None)

**Returns:** `BinarizationResult`

#### `trinarize_kmeans(vect, nstart=1, iter_max=10, na_rm=False, random_state=None)`

Trinarize a vector using k-Means clustering.

**Parameters:** Same as `binarize_kmeans`

**Returns:** `TrinarizationResult`

### Matrix Functions

#### `binarize_matrix(mat, method="kMeans", layer=None, **kwargs)`

Binarize a matrix row-by-row.

**Parameters:**
- `mat`: Input matrix (np.ndarray, pd.DataFrame, or AnnData)
- `method`: Discretization method (default: "kMeans")
- `layer`: For AnnData, which layer to use (default: None uses .X)
- `**kwargs`: Additional arguments for the binarization function

**Returns:** `pandas.DataFrame`

#### `trinarize_matrix(mat, method="kMeans", layer=None, **kwargs)`

Trinarize a matrix row-by-row.

**Parameters:** Same as `binarize_matrix`

**Returns:** `pandas.DataFrame`

## Examples

See the `examples/` directory for Jupyter notebooks demonstrating:
- Basic binarization and trinarization
- Working with gene expression matrices
- Integration with scanpy/AnnData workflows
- Visualization of results

## Relationship to R Package

This Python implementation provides the core k-Means-based methods from the original R BiTrinA package. The R package includes additional advanced methods (BASC, TASC) that may be added in future versions.

## Requirements

- Python >= 3.8
- NumPy >= 1.20.0
- Pandas >= 1.3.0
- scikit-learn >= 1.0.0
- anndata >= 0.8.0 (optional, for AnnData support)

## License

Artistic-2.0 (same as the original R package)

## Citation

If you use BiTrinA in your research, please cite the original R package:

```
Mundus S, Müssel C, Schmid F, Lausser L, Blätte TJ, Hopfensitz M, Kaftan P, Kestler HA (2025). 
BiTrinA: Binarization and Trinarization of One-Dimensional Data. 
R package version 1.3.2.
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Authors

Python implementation by the BiTrinA contributors.
Original R package by Stefan Mundus, Christoph Müssel, and others.
