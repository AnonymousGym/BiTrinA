"""
Unit tests for binarization functionality.
"""

import pytest
import numpy as np
from bitrina import binarize_kmeans
from bitrina.results import BinarizationResult


class TestBinarizeKMeans:
    """Tests for binarize_kmeans function."""
    
    def test_basic_binarization(self):
        """Test basic binarization with clear separation."""
        data = np.array([1.0, 1.1, 1.2, 5.0, 5.1, 5.2])
        result = binarize_kmeans(data, random_state=42)
        
        assert isinstance(result, BinarizationResult)
        assert len(result.binarized_measurements) == len(data)
        assert set(result.binarized_measurements) == {0, 1}
        assert result.method == "k-Means"
        assert 1.2 < result.threshold < 5.0
    
    def test_all_low_values_get_zero(self):
        """Test that lower cluster gets 0."""
        data = np.array([1.0, 1.1, 1.2, 5.0, 5.1, 5.2])
        result = binarize_kmeans(data, random_state=42)
        
        # First three values (low) should be 0
        assert np.all(result.binarized_measurements[:3] == 0)
        # Last three values (high) should be 1
        assert np.all(result.binarized_measurements[3:] == 1)
    
    def test_reproducibility_with_random_state(self):
        """Test that results are reproducible with random_state."""
        data = np.array([1.0, 1.5, 2.0, 5.0, 5.5, 6.0])
        result1 = binarize_kmeans(data, random_state=42)
        result2 = binarize_kmeans(data, random_state=42)
        
        np.testing.assert_array_equal(
            result1.binarized_measurements,
            result2.binarized_measurements
        )
        assert result1.threshold == result2.threshold
    
    def test_minimum_length_requirement(self):
        """Test that vector must have at least 3 entries."""
        with pytest.raises(ValueError, match="at least 3 entries"):
            binarize_kmeans([1.0, 2.0])
    
    def test_constant_vector_raises_error(self):
        """Test that constant vector raises error."""
        with pytest.raises(ValueError, match="constant"):
            binarize_kmeans([5.0, 5.0, 5.0, 5.0])
    
    def test_na_handling_with_na_rm_true(self):
        """Test NA removal when na_rm=True."""
        data = np.array([1.0, np.nan, 1.2, 5.0, np.nan, 5.2])
        result = binarize_kmeans(data, na_rm=True, random_state=42)
        
        # Should have 4 valid values
        assert len(result.original_measurements) == 4
        assert not np.any(np.isnan(result.original_measurements))
    
    def test_na_handling_with_na_rm_false(self):
        """Test that NA values raise error when na_rm=False."""
        data = np.array([1.0, np.nan, 1.2, 5.0, 5.2])
        with pytest.raises(ValueError, match="NA values"):
            binarize_kmeans(data, na_rm=False)
    
    def test_inf_values_raise_error(self):
        """Test that Inf values raise error."""
        data = np.array([1.0, 1.2, np.inf, 5.0, 5.2])
        with pytest.raises(ValueError, match="Inf values"):
            binarize_kmeans(data)
    
    def test_list_input(self):
        """Test that list input is handled correctly."""
        data = [1.0, 1.1, 1.2, 5.0, 5.1, 5.2]
        result = binarize_kmeans(data, random_state=42)
        
        assert isinstance(result, BinarizationResult)
        assert len(result.binarized_measurements) == len(data)
    
    def test_negative_nstart_raises_error(self):
        """Test that negative nstart raises error."""
        data = np.array([1.0, 1.1, 1.2, 5.0, 5.1, 5.2])
        with pytest.raises(ValueError, match="nstart"):
            binarize_kmeans(data, nstart=-1)
    
    def test_negative_iter_max_raises_error(self):
        """Test that negative iter_max raises error."""
        data = np.array([1.0, 1.1, 1.2, 5.0, 5.1, 5.2])
        with pytest.raises(ValueError, match="iter_max"):
            binarize_kmeans(data, iter_max=-1)
    
    def test_threshold_calculation(self):
        """Test that threshold is between clusters."""
        data = np.array([1.0, 1.0, 1.0, 10.0, 10.0, 10.0])
        result = binarize_kmeans(data, random_state=42)
        
        # Threshold should be between 1.0 and 10.0
        assert 1.0 < result.threshold < 10.0
        # Specifically, it should be close to midpoint
        assert 5.0 < result.threshold < 6.0
    
    def test_original_measurements_preserved(self):
        """Test that original measurements are preserved in result."""
        data = np.array([1.5, 2.0, 2.5, 7.0, 7.5, 8.0])
        result = binarize_kmeans(data, random_state=42)
        
        np.testing.assert_array_equal(result.original_measurements, data)
