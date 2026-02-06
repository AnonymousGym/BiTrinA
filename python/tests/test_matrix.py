"""
Unit tests for matrix operations.
"""

import pytest
import numpy as np
import pandas as pd
from bitrina import binarize_matrix, trinarize_matrix


class TestBinarizeMatrix:
    """Tests for binarize_matrix function."""
    
    def test_numpy_array_input(self):
        """Test binarization of numpy array."""
        mat = np.array([
            [1.0, 1.1, 5.0, 5.1],
            [2.0, 2.1, 8.0, 8.1]
        ])
        result = binarize_matrix(mat, random_state=42)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 6)  # 4 data cols + threshold + p_value
        assert 'threshold' in result.columns
        assert 'p_value' in result.columns
    
    def test_pandas_dataframe_input(self):
        """Test binarization of pandas DataFrame."""
        mat = pd.DataFrame({
            'Cell1': [1.0, 2.0],
            'Cell2': [1.1, 2.1],
            'Cell3': [5.0, 8.0],
            'Cell4': [5.1, 8.1]
        }, index=['Gene1', 'Gene2'])
        
        result = binarize_matrix(mat, random_state=42)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.index) == ['Gene1', 'Gene2']
        assert 'Cell1' in result.columns
        assert 'threshold' in result.columns
    
    def test_binarized_values_correct(self):
        """Test that binarized values are correct."""
        mat = np.array([
            [1.0, 1.0, 5.0, 5.0]
        ])
        result = binarize_matrix(mat, random_state=42)
        
        # First two should be 0, last two should be 1
        assert result.iloc[0, 0] == 0
        assert result.iloc[0, 1] == 0
        assert result.iloc[0, 2] == 1
        assert result.iloc[0, 3] == 1
    
    def test_threshold_values_calculated(self):
        """Test that threshold values are calculated for each row."""
        mat = np.array([
            [1.0, 1.0, 5.0, 5.0],
            [2.0, 2.0, 10.0, 10.0]
        ])
        result = binarize_matrix(mat, random_state=42)
        
        # Check thresholds are between clusters
        assert 1.0 < result.iloc[0]['threshold'] < 5.0
        assert 2.0 < result.iloc[1]['threshold'] < 10.0
    
    def test_unknown_method_raises_error(self):
        """Test that unknown method raises error."""
        mat = np.array([[1.0, 2.0, 3.0, 4.0]])
        with pytest.raises(ValueError, match="Unknown method"):
            binarize_matrix(mat, method="UnknownMethod")


class TestTrinarizeMatrix:
    """Tests for trinarize_matrix function."""
    
    def test_numpy_array_input(self):
        """Test trinarization of numpy array."""
        mat = np.array([
            [1.0, 1.1, 3.0, 3.1, 5.0, 5.1],
            [2.0, 2.1, 6.0, 6.1, 10.0, 10.1]
        ])
        result = trinarize_matrix(mat, random_state=42)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 9)  # 6 data cols + threshold1 + threshold2 + p_value
        assert 'threshold1' in result.columns
        assert 'threshold2' in result.columns
        assert 'p_value' in result.columns
    
    def test_pandas_dataframe_input(self):
        """Test trinarization of pandas DataFrame."""
        mat = pd.DataFrame({
            'Cell1': [1.0, 2.0],
            'Cell2': [1.1, 2.1],
            'Cell3': [3.0, 6.0],
            'Cell4': [3.1, 6.1],
            'Cell5': [5.0, 10.0],
            'Cell6': [5.1, 10.1]
        }, index=['Gene1', 'Gene2'])
        
        result = trinarize_matrix(mat, random_state=42)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.index) == ['Gene1', 'Gene2']
        assert 'Cell1' in result.columns
        assert 'threshold1' in result.columns
        assert 'threshold2' in result.columns
    
    def test_trinarized_values_correct(self):
        """Test that trinarized values are correct."""
        mat = np.array([
            [1.0, 1.0, 3.0, 3.0, 5.0, 5.0]
        ])
        result = trinarize_matrix(mat, random_state=42)
        
        # First two should be 0, middle two should be 1, last two should be 2
        assert result.iloc[0, 0] == 0
        assert result.iloc[0, 1] == 0
        assert result.iloc[0, 2] == 1
        assert result.iloc[0, 3] == 1
        assert result.iloc[0, 4] == 2
        assert result.iloc[0, 5] == 2
    
    def test_threshold_values_calculated(self):
        """Test that threshold values are calculated for each row."""
        mat = np.array([
            [1.0, 1.0, 5.0, 5.0, 10.0, 10.0]
        ])
        result = trinarize_matrix(mat, random_state=42)
        
        # Check thresholds are ordered and between clusters
        thresh1 = result.iloc[0]['threshold1']
        thresh2 = result.iloc[0]['threshold2']
        assert 1.0 < thresh1 < 5.0
        assert 5.0 < thresh2 < 10.0
        assert thresh1 < thresh2
    
    def test_unknown_method_raises_error(self):
        """Test that unknown method raises error."""
        mat = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
        with pytest.raises(ValueError, match="Unknown method"):
            trinarize_matrix(mat, method="UnknownMethod")
