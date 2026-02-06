"""
Basic usage examples for BiTrinA Python package.

This script demonstrates the core functionality of binarization and trinarization
using the k-Means method.
"""

import numpy as np
import pandas as pd
from bitrina import binarize_kmeans, trinarize_kmeans, binarize_matrix, trinarize_matrix


def example_1_basic_binarization():
    """Example 1: Basic binarization of a vector."""
    print("=" * 60)
    print("Example 1: Basic Binarization")
    print("=" * 60)
    
    # Create sample data with two clear groups
    data = np.array([1.2, 1.5, 1.3, 1.4, 5.6, 5.9, 5.4, 5.7])
    print(f"Original data: {data}")
    
    # Binarize using k-Means
    result = binarize_kmeans(data, random_state=42)
    print(f"Binarized:     {result.binarized_measurements}")
    print(f"Threshold:     {result.threshold:.4f}")
    print(f"Method:        {result.method}")
    print()


def example_2_basic_trinarization():
    """Example 2: Basic trinarization of a vector."""
    print("=" * 60)
    print("Example 2: Basic Trinarization")
    print("=" * 60)
    
    # Create sample data with three clear groups
    data = np.array([1.1, 1.2, 1.0, 3.5, 3.6, 3.4, 6.1, 6.2, 6.0])
    print(f"Original data: {data}")
    
    # Trinarize using k-Means
    result = trinarize_kmeans(data, random_state=42)
    print(f"Trinarized:    {result.trinarized_measurements}")
    print(f"Threshold 1:   {result.threshold1:.4f}")
    print(f"Threshold 2:   {result.threshold2:.4f}")
    print(f"Method:        {result.method}")
    print()


def example_3_gene_expression_matrix():
    """Example 3: Binarizing a gene expression matrix."""
    print("=" * 60)
    print("Example 3: Gene Expression Matrix Binarization")
    print("=" * 60)
    
    # Simulate gene expression data
    # 4 genes × 8 cells
    np.random.seed(42)
    expression_data = np.array([
        [1.2, 1.5, 1.3, 1.4, 5.6, 5.9, 5.4, 5.7],  # Gene A: low/high
        [0.1, 0.2, 0.15, 0.18, 2.1, 2.3, 2.0, 2.2],  # Gene B: off/on
        [3.1, 3.2, 3.0, 3.15, 7.5, 7.8, 7.6, 7.7],  # Gene C: medium/high
        [2.0, 2.1, 1.9, 2.05, 2.15, 2.2, 2.0, 2.1],  # Gene D: constant-ish
    ])
    
    gene_names = ["Gene_A", "Gene_B", "Gene_C", "Gene_D"]
    cell_names = [f"Cell_{i+1}" for i in range(8)]
    
    # Create DataFrame
    df = pd.DataFrame(expression_data, index=gene_names, columns=cell_names)
    print("Original expression data:")
    print(df)
    print()
    
    # Binarize the matrix
    binarized = binarize_matrix(df, random_state=42)
    print("Binarized expression data:")
    print(binarized[cell_names])  # Show only cell columns
    print()
    print("Thresholds for each gene:")
    print(binarized[["threshold"]])
    print()


def example_4_trinarize_matrix():
    """Example 4: Trinarizing a gene expression matrix."""
    print("=" * 60)
    print("Example 4: Gene Expression Matrix Trinarization")
    print("=" * 60)
    
    # Simulate gene expression with three levels
    np.random.seed(42)
    expression_data = np.array([
        [1.0, 1.1, 1.2, 3.5, 3.6, 3.7, 6.1, 6.2, 6.3],  # Gene A
        [0.5, 0.6, 0.55, 2.5, 2.6, 2.55, 5.5, 5.6, 5.55],  # Gene B
    ])
    
    gene_names = ["Gene_A", "Gene_B"]
    cell_names = [f"Cell_{i+1}" for i in range(9)]
    
    df = pd.DataFrame(expression_data, index=gene_names, columns=cell_names)
    print("Original expression data:")
    print(df)
    print()
    
    # Trinarize the matrix
    trinarized = trinarize_matrix(df, random_state=42)
    print("Trinarized expression data:")
    print(trinarized[cell_names])
    print()
    print("Thresholds for each gene:")
    print(trinarized[["threshold1", "threshold2"]])
    print()


def example_5_with_na_handling():
    """Example 5: Handling missing values."""
    print("=" * 60)
    print("Example 5: Handling Missing Values")
    print("=" * 60)
    
    # Create data with NaN values
    data = np.array([1.2, np.nan, 1.3, 5.6, 5.9, np.nan])
    print(f"Original data with NaN: {data}")
    
    # Remove NaN and binarize
    result = binarize_kmeans(data, na_rm=True, random_state=42)
    print(f"Binarized (na_rm=True): {result.binarized_measurements}")
    print(f"Threshold:              {result.threshold:.4f}")
    print(f"Number of valid values: {len(result.original_measurements)}")
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("BiTrinA Python Package - Usage Examples")
    print("=" * 60 + "\n")
    
    example_1_basic_binarization()
    example_2_basic_trinarization()
    example_3_gene_expression_matrix()
    example_4_trinarize_matrix()
    example_5_with_na_handling()
    
    print("=" * 60)
    print("Examples completed!")
    print("=" * 60)
