"""
Unit tests for trinarization functionality.
"""

import pytest
import numpy as np
from bitrina import trinarize_kmeans
from bitrina.results import TrinarizationResult


class TestTrinarizeKMeans:
    """Tests for trinarize_kmeans function."""
    
    def test_basic_trinarization(self):
        """Test basic trinarization with clear separation."""
        data = np.array([1.0, 1.1, 3.0, 3.1, 5.0, 5.1])
        result = trinarize_kmeans(data, random_state=42)
        
        assert isinstance(result, TrinarizationResult)
        assert len(result.trinarized_measurements) == len(data)
        assert set(result.trinarized_measurements) == {0, 1, 2}
        assert result.method == "k-Means"
        assert result.threshold1 < result.threshold2
    
    def test_correct_ordering(self):
        """Test that clusters are ordered correctly (0 < 1 < 2)."""
        data = np.array([1.0, 1.1, 3.0, 3.1, 5.0, 5.1])
        result = trinarize_kmeans(data, random_state=42)
        
        # Lowest values should be 0
        assert np.all(result.trinarized_measurements[:2] == 0)
        # Middle values should be 1
        assert np.all(result.trinarized_measurements[2:4] == 1)
        # Highest values should be 2
        assert np.all(result.trinarized_measurements[4:] == 2)
    
    def test_reproducibility_with_random_state(self):
        """Test that results are reproducible with random_state."""
        data = np.array([1.0, 1.5, 3.0, 3.5, 5.0, 5.5])
        result1 = trinarize_kmeans(data, random_state=42)
        result2 = trinarize_kmeans(data, random_state=42)
        
        np.testing.assert_array_equal(
            result1.trinarized_measurements,
            result2.trinarized_measurements
        )
        assert result1.threshold1 == result2.threshold1
        assert result1.threshold2 == result2.threshold2
    
    def test_minimum_length_requirement(self):
        """Test that vector must have at least 3 entries."""
        with pytest.raises(ValueError, match="at least 3 entries"):
            trinarize_kmeans([1.0, 2.0])
    
    def test_constant_vector_raises_error(self):
        """Test that constant vector raises error."""
        with pytest.raises(ValueError, match="constant"):
            trinarize_kmeans([5.0, 5.0, 5.0, 5.0])
    
    def test_na_handling_with_na_rm_true(self):
        """Test NA removal when na_rm=True."""
        data = np.array([1.0, np.nan, 3.0, np.nan, 5.0, 5.1])
        result = trinarize_kmeans(data, na_rm=True, random_state=42)
        
        # Should have 4 valid values
        assert len(result.original_measurements) == 4
        assert not np.any(np.isnan(result.original_measurements))
    
    def test_na_handling_with_na_rm_false(self):
        """Test that NA values raise error when na_rm=False."""
        data = np.array([1.0, np.nan, 3.0, 5.0])
        with pytest.raises(ValueError, match="NA values"):
            trinarize_kmeans(data, na_rm=False)
    
    def test_inf_values_raise_error(self):
        """Test that Inf values raise error."""
        data = np.array([1.0, 3.0, np.inf, 5.0])
        with pytest.raises(ValueError, match="Inf values"):
            trinarize_kmeans(data)
    
    def test_list_input(self):
        """Test that list input is handled correctly."""
        data = [1.0, 1.1, 3.0, 3.1, 5.0, 5.1]
        result = trinarize_kmeans(data, random_state=42)
        
        assert isinstance(result, TrinarizationResult)
        assert len(result.trinarized_measurements) == len(data)
    
    def test_threshold_ordering(self):
        """Test that thresholds are ordered correctly."""
        data = np.array([1.0, 1.0, 5.0, 5.0, 10.0, 10.0])
        result = trinarize_kmeans(data, random_state=42)
        
        # threshold1 should separate low and mid
        assert 1.0 < result.threshold1 < 5.0
        # threshold2 should separate mid and high
        assert 5.0 < result.threshold2 < 10.0
        # threshold1 < threshold2
        assert result.threshold1 < result.threshold2
    
    def test_original_measurements_preserved(self):
        """Test that original measurements are preserved in result."""
        data = np.array([1.5, 2.0, 4.0, 4.5, 7.0, 7.5])
        result = trinarize_kmeans(data, random_state=42)
        
        np.testing.assert_array_equal(result.original_measurements, data)
    
    def test_negative_nstart_raises_error(self):
        """Test that negative nstart raises error."""
        data = np.array([1.0, 3.0, 5.0])
        with pytest.raises(ValueError, match="nstart"):
            trinarize_kmeans(data, nstart=-1)
    
    def test_negative_iter_max_raises_error(self):
        """Test that negative iter_max raises error."""
        data = np.array([1.0, 3.0, 5.0])
        with pytest.raises(ValueError, match="iter_max"):
            trinarize_kmeans(data, iter_max=-1)
