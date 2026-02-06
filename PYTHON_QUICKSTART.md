# BiTrinA Python Implementation - Quick Start

## What is BiTrinA?

BiTrinA (Binarization and Trinarization) converts continuous numerical data into discrete categories:
- **Binarization**: Two states (0 and 1) - e.g., "off" and "on"
- **Trinarization**: Three states (0, 1, and 2) - e.g., "low", "medium", and "high"

This Python implementation is specifically designed for **single-cell RNA-seq data** using the AnnData format.

## Installation

```bash
cd python
pip install -e .

# For single-cell RNA-seq support
pip install -e ".[anndata]"
```

## Quick Example: Single-Cell RNA-seq

```python
import scanpy as sc
from bitrina import binarize_matrix

# Load your scRNA-seq data
adata = sc.read_h5ad("your_data.h5ad")

# Preprocess (optional)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Binarize gene expression
# Note: Use .T to transpose (process genes independently across cells)
binarized_df = binarize_matrix(adata.T, method="kMeans")

# Store results back in AnnData
adata.layers['binarized'] = binarized_df.iloc[:, :-2].T.values

# Now you can use binarized expression for downstream analysis
```

## Why Transpose `adata.T`?

- Standard AnnData shape: `(n_cells, n_genes)`
- BiTrinA processes each **row** independently
- We want to process each **gene** independently across all cells
- Solution: Transpose with `adata.T` to get `(n_genes, n_cells)` shape

## Input and Output

### Input Formats
1. **AnnData objects** (for scRNA-seq) - must transpose with `.T`
2. **NumPy arrays** - shape `(n_features, n_samples)`
3. **Pandas DataFrames** - rows = features, columns = samples

### Output Format
- **pandas DataFrame** with:
  - Discretized values for each sample (column)
  - Metadata columns: `threshold` (or `threshold1`, `threshold2`), `p_value`
  - Preserves row and column names from input

## More Examples

### Basic Binarization of a Vector

```python
from bitrina import binarize_kmeans
import numpy as np

data = np.array([1.2, 1.5, 1.3, 5.6, 5.9, 5.4])
result = binarize_kmeans(data, random_state=42)

print(result.binarized_measurements)  # [0 0 0 1 1 1]
print(result.threshold)  # 3.45
```

### Trinarization for Three Expression Levels

```python
from bitrina import trinarize_matrix

# Trinarize to detect: off (0), low (1), high (2)
trinarized_df = trinarize_matrix(adata.T, method="kMeans")

# Store in AnnData
adata.layers['trinarized'] = trinarized_df.iloc[:, :-3].T.values
```

### Working with Gene Expression Matrix

```python
from bitrina import binarize_matrix
import pandas as pd

# Gene expression: genes × cells
expression = pd.DataFrame({
    'Cell1': [1.2, 0.1, 3.1],
    'Cell2': [1.5, 0.2, 3.2],
    'Cell3': [5.6, 2.1, 7.5],
    'Cell4': [5.9, 2.3, 7.8]
}, index=['GeneA', 'GeneB', 'GeneC'])

binarized = binarize_matrix(expression, random_state=42)
print(binarized)
```

## Documentation

- **Full Documentation**: `python/README.md`
- **Input/Output Specs**: `python/INPUT_OUTPUT_SPEC.md`
- **Implementation Summary**: `python/SUMMARY.md`
- **Examples**: `python/examples/`

## Run Examples

```bash
cd python

# Basic usage examples
python examples/basic_usage.py

# AnnData integration examples
python examples/anndata_integration.py
```

## Run Tests

```bash
cd python
pytest tests/ -v
```

All 36 tests should pass ✓

## Common Use Cases

### 1. Identify "On/Off" Genes

```python
binarized = binarize_matrix(adata.T)
on_fraction = (binarized.iloc[:, :-2] == 1).mean(axis=1)
always_on_genes = on_fraction[on_fraction > 0.8].index
```

### 2. Find Bimodal Genes

```python
binarized = binarize_matrix(adata.T)
on_fraction = (binarized.iloc[:, :-2] == 1).mean(axis=1)
bimodal = on_fraction[(on_fraction > 0.3) & (on_fraction < 0.7)]
```

### 3. Cell Type Markers

```python
for cell_type in adata.obs['cell_type'].unique():
    cells = adata.obs['cell_type'] == cell_type
    type_expression = binarized.loc[:, cells].mean(axis=1)
    markers = type_expression[type_expression > 0.8].index
    print(f"{cell_type}: {len(markers)} markers")
```

## API Reference

### Vector Functions

- `binarize_kmeans(vect, nstart=1, iter_max=10, na_rm=False, random_state=None)`
- `trinarize_kmeans(vect, nstart=1, iter_max=10, na_rm=False, random_state=None)`

### Matrix Functions

- `binarize_matrix(mat, method="kMeans", layer=None, **kwargs)`
- `trinarize_matrix(mat, method="kMeans", layer=None, **kwargs)`

## Methods Available

Currently implemented:
- **k-Means** (`method="kMeans"`): Fast clustering-based approach

Future methods (from R package):
- BASC A/B: Advanced scale-space binarization
- TASC A/B: Advanced scale-space trinarization

## Requirements

- Python >= 3.8
- numpy >= 1.20.0
- pandas >= 1.3.0
- scikit-learn >= 1.0.0
- anndata >= 0.8.0 (optional, for scRNA-seq)

## Getting Help

1. Check the documentation in `python/README.md`
2. Look at examples in `python/examples/`
3. Read the detailed I/O spec in `python/INPUT_OUTPUT_SPEC.md`
4. Review the implementation summary in `python/SUMMARY.md`

## Citation

If you use BiTrinA in your research, please cite:

```
Mundus S, Müssel C, Schmid F, Lausser L, Blätte TJ, Hopfensitz M, 
Kaftan P, Kestler HA (2025). BiTrinA: Binarization and Trinarization 
of One-Dimensional Data. R package version 1.3.2.
```

## License

Artistic-2.0 (same as original R package)
