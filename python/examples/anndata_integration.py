"""
Example demonstrating integration with AnnData for single-cell RNA-seq analysis.

Note: This example requires anndata to be installed:
    pip install anndata
"""

import numpy as np


def create_example_adata():
    """Create a simple example AnnData object for demonstration."""
    try:
        import anndata as ad
    except ImportError:
        print("This example requires anndata. Install with: pip install anndata")
        return None
    
    # Create synthetic single-cell data
    # 100 cells × 50 genes
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    # Simulate expression with some structure
    # - 20 genes highly expressed in first 50 cells
    # - 15 genes highly expressed in last 50 cells
    # - 15 genes with moderate expression everywhere
    
    expression = np.random.lognormal(0, 0.5, size=(n_cells, n_genes))
    
    # Add differential expression patterns
    expression[:50, :20] *= 5  # High in first group
    expression[50:, 20:35] *= 5  # High in second group
    
    # Create AnnData object
    adata = ad.AnnData(X=expression)
    
    # Add cell metadata
    adata.obs['cell_type'] = ['Type_A'] * 50 + ['Type_B'] * 50
    adata.obs['batch'] = np.random.choice(['Batch1', 'Batch2'], n_cells)
    
    # Add gene metadata
    adata.var['gene_name'] = [f'Gene_{i+1}' for i in range(n_genes)]
    adata.var['highly_variable'] = [True] * 35 + [False] * 15
    
    return adata


def example_anndata_binarization():
    """Example: Binarize gene expression in AnnData object."""
    print("=" * 60)
    print("Example: AnnData Binarization")
    print("=" * 60)
    
    # Create example data
    adata = create_example_adata()
    if adata is None:
        return
    
    print(f"AnnData shape: {adata.shape} (cells × genes)")
    print(f"Cell types: {adata.obs['cell_type'].unique()}")
    print()
    
    # Import BiTrinA
    from bitrina import binarize_matrix
    
    # Binarize gene expression
    # Note: We transpose because binarize_matrix processes rows independently
    # and we want to process each gene independently across cells
    print("Binarizing gene expression...")
    binarized_df = binarize_matrix(adata.T, method="kMeans", random_state=42)
    
    print(f"Binarized data shape: {binarized_df.shape}")
    print(f"Columns: {list(binarized_df.columns[:5])} ... (first 5 cells shown)")
    print()
    
    # Show results for first few genes
    print("First 5 genes - binarized expression for first 8 cells:")
    print(binarized_df.iloc[:5, :8])
    print()
    
    print("Thresholds for first 5 genes:")
    print(binarized_df.iloc[:5, -2:])  # Last 2 columns are threshold and p_value
    print()
    
    # Store back in AnnData as a layer
    # Remove threshold and p_value columns first
    binarized_values = binarized_df.iloc[:, :-2].T.values
    adata.layers['binarized'] = binarized_values
    
    print(f"Stored binarized data in adata.layers['binarized']")
    print(f"Layer shape: {adata.layers['binarized'].shape}")
    print()


def example_anndata_trinarization():
    """Example: Trinarize gene expression in AnnData object."""
    print("=" * 60)
    print("Example: AnnData Trinarization")
    print("=" * 60)
    
    # Create example data
    adata = create_example_adata()
    if adata is None:
        return
    
    # Import BiTrinA
    from bitrina import trinarize_matrix
    
    # Trinarize gene expression
    print("Trinarizing gene expression...")
    trinarized_df = trinarize_matrix(adata.T, method="kMeans", random_state=42)
    
    print(f"Trinarized data shape: {trinarized_df.shape}")
    print()
    
    # Show results for first few genes
    print("First 5 genes - trinarized expression for first 8 cells:")
    print(trinarized_df.iloc[:5, :8])
    print()
    
    print("Thresholds for first 5 genes:")
    print(trinarized_df.iloc[:5, -3:])  # Last 3 columns are threshold1, threshold2, p_value
    print()
    
    # Store back in AnnData
    trinarized_values = trinarized_df.iloc[:, :-3].T.values
    adata.layers['trinarized'] = trinarized_values
    
    print(f"Stored trinarized data in adata.layers['trinarized']")
    print()


def example_specific_layer():
    """Example: Working with specific layers in AnnData."""
    print("=" * 60)
    print("Example: Working with AnnData Layers")
    print("=" * 60)
    
    # Create example data
    adata = create_example_adata()
    if adata is None:
        return
    
    # Add a normalized layer
    import numpy as np
    adata.layers['normalized'] = np.log1p(adata.X)
    
    print(f"Available layers: {list(adata.layers.keys())}")
    print()
    
    # Import BiTrinA
    from bitrina import binarize_matrix
    
    # Binarize the normalized layer
    print("Binarizing normalized layer...")
    binarized_df = binarize_matrix(
        adata.T,
        method="kMeans",
        layer="normalized",
        random_state=42
    )
    
    print(f"Binarized normalized data shape: {binarized_df.shape}")
    print("First 3 genes:")
    print(binarized_df.iloc[:3, :8])
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("BiTrinA + AnnData Integration Examples")
    print("=" * 60 + "\n")
    
    example_anndata_binarization()
    example_anndata_trinarization()
    example_specific_layer()
    
    print("=" * 60)
    print("Examples completed!")
    print("=" * 60)
