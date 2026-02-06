# BiTrinA Python Implementation - Summary

## Overview

This document provides a comprehensive summary of the Python implementation of BiTrinA (Binarization and Trinarization), created to work with single-cell RNA-seq data in the AnnData format.

## What Was Implemented

### 1. Core Package Structure (`python/bitrina/`)

A complete Python package with the following modules:

#### **`results.py`** - Result Data Classes
- `BinarizationResult`: Stores binarization outputs
- `TrinarizationResult`: Stores trinarization outputs
- Both include: original data, discretized values, thresholds, p-values, and method names

#### **`binarize.py`** - Binarization Functions
- `binarize_kmeans()`: Convert continuous data to binary (0/1) using k-Means clustering
- Handles vectors (1D arrays)
- Comprehensive input validation and error handling
- Support for NA value handling

#### **`trinarize.py`** - Trinarization Functions
- `trinarize_kmeans()`: Convert continuous data to ternary (0/1/2) using k-Means clustering
- Handles vectors (1D arrays)
- Automatically orders clusters from low to high
- Support for NA value handling

#### **`matrix.py`** - Batch Processing for Matrices
- `binarize_matrix()`: Process entire matrices row-by-row
- `trinarize_matrix()`: Process entire matrices row-by-row
- **Full AnnData support**: Handles single-cell RNA-seq data
- Support for multiple input formats:
  - NumPy arrays
  - Pandas DataFrames
  - **AnnData objects** (with automatic transpose handling)
- Layer support for AnnData
- Returns pandas DataFrames with metadata columns

### 2. Documentation

#### **`README.md`** - User Guide
- Installation instructions
- Quick start examples
- API reference
- Usage patterns for different data types
- Examples for scRNA-seq workflows

#### **`INPUT_OUTPUT_SPEC.md`** - Detailed I/O Specifications
- Comprehensive guide for working with scRNA-seq data
- Explains AnnData format and transpose requirements
- Common use cases and workflows
- Performance considerations
- Integration with Scanpy

### 3. Examples (`python/examples/`)

#### **`basic_usage.py`** - Core Functionality Examples
- Example 1: Basic binarization of vectors
- Example 2: Basic trinarization of vectors
- Example 3: Gene expression matrix binarization
- Example 4: Gene expression matrix trinarization
- Example 5: Handling missing values

#### **`anndata_integration.py`** - scRNA-seq Examples
- Creating example AnnData objects
- Binarizing gene expression across cells
- Trinarizing gene expression
- Working with AnnData layers
- Storing results back in AnnData

### 4. Tests (`python/tests/`)

Comprehensive test suite with **36 tests** covering:

#### **`test_binarize.py`** (13 tests)
- Basic functionality
- Reproducibility
- Input validation
- NA/Inf value handling
- Threshold calculations

#### **`test_trinarize.py`** (13 tests)
- Basic functionality
- Cluster ordering
- Reproducibility
- Input validation
- Threshold calculations

#### **`test_matrix.py`** (10 tests)
- Multiple input formats
- Column/row name preservation
- Value correctness
- Threshold calculations

**All tests pass ✓**

## Input and Output Formats

### Input

BiTrinA accepts three main input formats:

1. **AnnData objects** (recommended for scRNA-seq)
   ```python
   binarize_matrix(adata.T, method="kMeans")
   ```
   - Must transpose (`adata.T`) to process genes independently
   - Can specify layers: `layer="log1p"`
   
2. **NumPy arrays**
   ```python
   binarize_matrix(np.array([[1,2,3], [4,5,6]]))
   ```
   - Shape: (n_features, n_samples)
   - Rows processed independently

3. **Pandas DataFrames**
   ```python
   binarize_matrix(df)
   ```
   - Preserves row/column names
   - Index = features, Columns = samples

### Output

Matrix operations return **pandas DataFrames**:

**Binarization output:**
- Data columns: Binarized values (0 or 1)
- Metadata columns: `threshold`, `p_value`

**Trinarization output:**
- Data columns: Trinarized values (0, 1, or 2)
- Metadata columns: `threshold1`, `threshold2`, `p_value`

Vector operations return **Result objects**:
- `BinarizationResult` or `TrinarizationResult`
- Attributes: original_measurements, binarized/trinarized_measurements, thresholds, p_value, method

## Key Features

### 1. scRNA-seq Integration

**Full AnnData Support:**
- Automatic handling of transposed AnnData
- Layer support for different normalization states
- Preserves cell and gene names
- Easy storage back into AnnData layers

**Example Workflow:**
```python
import scanpy as sc
from bitrina import binarize_matrix

# Load data
adata = sc.read_h5ad("pbmc3k.h5ad")

# Preprocess
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Binarize (note the transpose)
binarized_df = binarize_matrix(adata.T, method="kMeans")

# Store results
adata.layers['binarized'] = binarized_df.iloc[:, :-2].T.values
```

### 2. Robust Error Handling

- Input validation (type, size, values)
- Graceful handling of constant vectors
- NA/NaN value support with `na_rm` parameter
- Clear error messages

### 3. Method Implementation

Currently implements **k-Means clustering**:
- Simple and fast
- No external C dependencies (unlike R version's BASC/TASC)
- Uses scikit-learn's optimized implementation
- Reproducible with `random_state` parameter

### 4. Flexibility

- Works with vectors, matrices, and DataFrames
- Batch processing of multiple features
- Optional parameters for customization
- Extensible design for future methods

## Installation

```bash
cd python
pip install -e .

# For scRNA-seq support
pip install -e ".[anndata]"

# For development
pip install -e ".[dev]"
```

## Dependencies

**Core:**
- numpy >= 1.20.0
- pandas >= 1.3.0
- scikit-learn >= 1.0.0

**Optional:**
- anndata >= 0.8.0 (for scRNA-seq support)

**Development:**
- pytest >= 7.0.0
- pytest-cov >= 3.0.0

## Usage Examples

### Basic Binarization

```python
from bitrina import binarize_kmeans

data = [1.2, 1.5, 1.3, 5.6, 5.9, 5.4]
result = binarize_kmeans(data, random_state=42)
print(result.binarized_measurements)  # [0 0 0 1 1 1]
print(result.threshold)  # 3.45
```

### Gene Expression Matrix

```python
from bitrina import binarize_matrix
import pandas as pd

# Gene expression: 3 genes × 6 cells
expression = pd.DataFrame({
    'Cell1': [1.2, 0.1, 3.1],
    'Cell2': [1.5, 0.2, 3.2],
    'Cell3': [5.6, 2.1, 7.5],
    'Cell4': [5.9, 2.3, 7.8]
}, index=['GeneA', 'GeneB', 'GeneC'])

binarized = binarize_matrix(expression)
print(binarized)
```

### Single-Cell RNA-seq

```python
import scanpy as sc
from bitrina import binarize_matrix

adata = sc.read_h5ad("data.h5ad")
sc.pp.log1p(adata)

# Binarize all genes across cells
binarized_df = binarize_matrix(adata.T, method="kMeans")

# Store in AnnData
adata.layers['binarized'] = binarized_df.iloc[:, :-2].T.values
```

## Comparison with R Package

### Similarities
- Same conceptual approach (discretization of continuous data)
- k-Means method implemented in both
- Similar API design
- Comprehensive result objects

### Differences

| Feature | R Package | Python Package |
|---------|-----------|----------------|
| **Methods** | k-Means, BASC A/B, TASC A/B | k-Means only (initially) |
| **Dependencies** | Requires C compilation | Pure Python (scikit-learn) |
| **Input Format** | Vectors, matrices | Vectors, matrices, DataFrames, AnnData |
| **Primary Use Case** | General-purpose | Optimized for scRNA-seq |
| **Integration** | Base R, Bioconductor | Scanpy, AnnData ecosystem |

## Testing

Run tests:
```bash
cd python
pytest tests/ -v
```

Run with coverage:
```bash
pytest tests/ --cov=bitrina --cov-report=html
```

## Future Extensions

Potential additions:
1. **BASC methods**: Advanced scale-space based binarization
2. **TASC methods**: Advanced scale-space based trinarization
3. **Statistical tests**: Implement dip test for p-values
4. **Visualization**: Plot functions for showing discretization
5. **Batch optimization**: Parallel processing for large datasets
6. **Sparse matrix support**: Direct handling of sparse AnnData matrices

## File Structure

```
python/
├── bitrina/
│   ├── __init__.py          # Package initialization
│   ├── binarize.py          # Binarization methods
│   ├── trinarize.py         # Trinarization methods
│   ├── matrix.py            # Matrix operations
│   └── results.py           # Result classes
├── examples/
│   ├── basic_usage.py       # Basic examples
│   └── anndata_integration.py  # scRNA-seq examples
├── tests/
│   ├── test_binarize.py     # Binarization tests
│   ├── test_trinarize.py    # Trinarization tests
│   └── test_matrix.py       # Matrix operation tests
├── README.md                # User guide
├── INPUT_OUTPUT_SPEC.md     # I/O specifications
├── setup.py                 # Package configuration
└── requirements.txt         # Dependencies
```

## License

Artistic-2.0 (same as original R package)

## Citation

If you use BiTrinA in your research, please cite the original R package:

```
Mundus S, Müssel C, Schmid F, Lausser L, Blätte TJ, Hopfensitz M, 
Kaftan P, Kestler HA (2025). BiTrinA: Binarization and Trinarization 
of One-Dimensional Data. R package version 1.3.2.
```

## Contributors

Python implementation: BiTrinA Contributors
Original R package: Stefan Mundus, Christoph Müssel, Hans A. Kestler, and others

---

**Implementation Date:** February 2026
**Version:** 0.1.0
**Status:** Production-ready for k-Means methods
