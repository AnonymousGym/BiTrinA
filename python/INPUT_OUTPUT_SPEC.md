# Input and Output Specifications for scRNA-seq Data

## Overview

BiTrinA's Python implementation is designed to work seamlessly with single-cell RNA-seq (scRNA-seq) data, particularly data stored in the **AnnData** format used by Scanpy and other popular single-cell analysis tools.

## Input Formats

### 1. AnnData Objects (Recommended for scRNA-seq)

**Structure:**
- AnnData objects store scRNA-seq data in a structured format
- Shape: `(n_obs, n_vars)` where `n_obs` = cells and `n_vars` = genes
- Main data matrix: `adata.X`
- Alternative data layers: `adata.layers['layer_name']`

**Example:**
```python
import scanpy as sc
from bitrina import binarize_matrix

# Load scRNA-seq data
adata = sc.read_h5ad("pbmc3k.h5ad")
# adata.shape: (2700 cells, 32738 genes)

# Binarize gene expression across cells
# Transpose to get genes as rows (one row per gene)
binarized_df = binarize_matrix(adata.T, method="kMeans")
```

**Key Points:**
- **Must transpose** (`adata.T`) because BiTrinA processes each row independently
- After transposition: shape becomes `(n_genes, n_cells)`
- Each gene is binarized/trinarized independently across all cells

### 2. NumPy Arrays

**Structure:**
- 2D array with shape `(n_features, n_samples)`
- Rows = features (e.g., genes)
- Columns = samples (e.g., cells)

**Example:**
```python
import numpy as np
from bitrina import binarize_matrix

# Gene expression matrix: 100 genes × 500 cells
expression_data = np.random.lognormal(0, 1, size=(100, 500))
binarized_df = binarize_matrix(expression_data, method="kMeans")
```

### 3. Pandas DataFrames

**Structure:**
- Rows = features (genes), with optional index
- Columns = samples (cells), with optional column names

**Example:**
```python
import pandas as pd
from bitrina import binarize_matrix

# Gene expression DataFrame
df = pd.DataFrame(
    expression_data,
    index=[f"Gene_{i}" for i in range(100)],
    columns=[f"Cell_{i}" for i in range(500)]
)
binarized_df = binarize_matrix(df, method="kMeans")
```

## Output Format

### Matrix Operations Return pandas DataFrames

**Structure:**
- Rows: Features (genes) - preserves input row names/index
- Columns: 
  - Sample columns (cells) - preserves input column names
  - `threshold` or `threshold1`, `threshold2`: computed thresholds
  - `p_value`: statistical test result (currently None for k-Means)

**Example Output:**
```python
>>> binarized_df = binarize_matrix(adata.T, method="kMeans")
>>> binarized_df.shape
(32738, 2702)  # 32738 genes × (2700 cells + threshold + p_value)

>>> binarized_df.head()
          Cell_1  Cell_2  Cell_3  ...  threshold  p_value
Gene_1         0       0       1  ...      2.456     None
Gene_2         1       1       1  ...      0.123     None
Gene_3         0       1       0  ...      1.789     None
```

## Typical scRNA-seq Workflow

### Standard Workflow

```python
import scanpy as sc
from bitrina import binarize_matrix, trinarize_matrix

# 1. Load data
adata = sc.read_h5ad("your_data.h5ad")

# 2. Preprocess (optional but recommended)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 3. Binarize expression
# Note: Transpose to get genes as rows
binarized_df = binarize_matrix(adata.T, method="kMeans")

# 4. Extract binarized values (excluding threshold and p_value columns)
binarized_values = binarized_df.iloc[:, :-2]

# 5. Store back in AnnData (optional)
# Transpose back to (cells × genes) format
adata.layers['binarized'] = binarized_values.T.values

# 6. Use binarized data for downstream analysis
# e.g., clustering, differential expression, etc.
```

### Working with Specific Layers

```python
# If you have multiple processed versions in layers
adata.layers['log1p'] = np.log1p(adata.X)
adata.layers['scaled'] = sc.pp.scale(adata.X, copy=True)

# Binarize a specific layer
binarized_log = binarize_matrix(
    adata.T,
    method="kMeans",
    layer="log1p"
)
```

### Trinarization for Three Expression Levels

```python
# Trinarize to detect: off (0), low (1), high (2)
trinarized_df = trinarize_matrix(adata.T, method="kMeans")

# Store in AnnData
trinarized_values = trinarized_df.iloc[:, :-3]  # Remove threshold1, threshold2, p_value
adata.layers['trinarized'] = trinarized_values.T.values
```

## Data Organization Conventions

### Why Transpose?

BiTrinA processes each **row** independently. For scRNA-seq:
- Standard AnnData shape: `(cells, genes)`
- We want to process each gene independently across cells
- Solution: Use `adata.T` to get `(genes, cells)` shape

### Storage in AnnData Layers

```python
# After binarization/trinarization
adata.layers['binarized'] = binarized_values.T.values
adata.layers['trinarized'] = trinarized_values.T.values

# Access later
binarized_expression = adata.layers['binarized']  # Shape: (cells, genes)
```

## Input Requirements and Constraints

### Vector Operations (`binarize_kmeans`, `trinarize_kmeans`)

**Required:**
- Numeric values (int or float)
- At least 3 values
- Non-constant (must have variation)

**Optional:**
- Handle NA/NaN values with `na_rm=True`
- Cannot handle Inf values (will raise error)

### Matrix Operations (`binarize_matrix`, `trinarize_matrix`)

**Required:**
- 2D structure (matrix, DataFrame, or AnnData)
- Rows processed independently
- Each row must meet vector requirements

**Handling Failed Rows:**
- If a row fails (e.g., constant values), it's filled with NaN
- Processing continues for other rows

## Common Use Cases

### 1. Boolean Gene Expression Analysis

```python
# Determine if genes are "on" or "off"
binarized = binarize_matrix(adata.T, method="kMeans")
gene_on_fraction = (binarized.iloc[:, :-2] == 1).mean(axis=1)
print(gene_on_fraction.sort_values(ascending=False).head(10))
```

### 2. Expression Level Categories

```python
# Categorize into low/medium/high
trinarized = trinarize_matrix(adata.T, method="kMeans")
# Count cells in each category per gene
for level in [0, 1, 2]:
    adata.var[f'n_cells_level_{level}'] = (
        trinarized.iloc[:, :-3] == level
    ).sum(axis=1).values
```

### 3. Feature Selection

```python
# Find genes with bimodal expression
binarized = binarize_matrix(adata.T, method="kMeans")
# Genes where both states are common (30-70% cells in each)
on_fraction = (binarized.iloc[:, :-2] == 1).mean(axis=1)
bimodal_genes = on_fraction[(on_fraction > 0.3) & (on_fraction < 0.7)].index
```

### 4. Cell Type Marker Identification

```python
# Find genes with clear on/off patterns per cell type
binarized = binarize_matrix(adata.T, method="kMeans")

for cell_type in adata.obs['cell_type'].unique():
    cells = adata.obs['cell_type'] == cell_type
    type_specific = binarized.loc[:, cells].mean(axis=1)
    # Genes highly expressed in this cell type
    markers = type_specific[type_specific > 0.8].index
    print(f"{cell_type}: {len(markers)} marker genes")
```

## Performance Considerations

### Large Datasets

- Processing is row-wise, so time scales with number of genes
- For 30,000 genes: ~10-30 seconds on typical hardware
- Memory: Primarily determined by input data size
- Consider subsetting to highly variable genes if needed

### Optimization Tips

```python
# 1. Subset to highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata_subset = adata[:, adata.var['highly_variable']]
binarized = binarize_matrix(adata_subset.T, method="kMeans")

# 2. Process in batches for very large datasets
batch_size = 1000
results = []
for i in range(0, adata.n_vars, batch_size):
    batch = adata[:, i:i+batch_size]
    result = binarize_matrix(batch.T, method="kMeans")
    results.append(result)
```

## References

- AnnData documentation: https://anndata.readthedocs.io/
- Scanpy documentation: https://scanpy.readthedocs.io/
- Original BiTrinA R package: https://cran.r-project.org/package=BiTrinA
