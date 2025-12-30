"""Tests for BioactivityExtractor."""

import numpy as np
import pandas as pd
from chembl_pdb_linker.extractors.bioactivity import BioactivityExtractor


class TestBioactivityExtractor:
    """Test cases for BioactivityExtractor."""

    def test_filter_activities_by_type(self, config):
        """Test filtering activities by type."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_type": ["IC50", "Ki", "EC50", "Kd", "MIC", "GI50"],
                "standard_value": [100, 200, 300, 400, 500, 600],
                "standard_units": ["nM"] * 6,
                "confidence_score": [9] * 6,
            }
        )

        result = extractor.filter_activities(
            activities,
            activity_types=["IC50", "Ki"],
        )

        assert len(result) == 2
        assert set(result["standard_type"]) == {"IC50", "Ki"}

    def test_filter_activities_by_confidence(self, config):
        """Test filtering activities by confidence score."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_type": ["IC50"] * 4,
                "standard_value": [100, 200, 300, 400],
                "standard_units": ["nM"] * 4,
                "confidence_score": [9, 8, 7, 6],
            }
        )

        result = extractor.filter_activities(activities, min_confidence=8)

        assert len(result) == 2
        assert all(result["confidence_score"] >= 8)

    def test_filter_activities_removes_null_values(self, config):
        """Test that null activity values are removed."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_type": ["IC50"] * 3,
                "standard_value": [100, None, 300],
                "standard_units": ["nM"] * 3,
                "confidence_score": [9] * 3,
            }
        )

        result = extractor.filter_activities(activities)

        assert len(result) == 2
        assert result["standard_value"].notna().all()

    def test_standardize_units_nm_to_nm(self, config):
        """Test that nM values remain unchanged."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_value": [100.0, 200.0],
                "standard_units": ["nM", "nM"],
            }
        )

        result = extractor.standardize_units(activities, target_unit="nM")

        assert "standardized_value" in result.columns
        assert result["standardized_value"].tolist() == [100.0, 200.0]

    def test_standardize_units_um_to_nm(self, config):
        """Test conversion from uM to nM (multiply by 1000)."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_value": [1.0, 10.0],
                "standard_units": ["uM", "uM"],
            }
        )

        result = extractor.standardize_units(activities, target_unit="nM")

        # 1 uM = 1000 nM
        assert result["standardized_value"].tolist() == [1000.0, 10000.0]

    def test_standardize_units_mm_to_nm(self, config):
        """Test conversion from mM to nM (multiply by 1,000,000)."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_value": [0.001],
                "standard_units": ["mM"],
            }
        )

        result = extractor.standardize_units(activities, target_unit="nM")

        # 0.001 mM = 1000 nM
        assert result["standardized_value"].tolist() == [1000.0]

    def test_standardize_units_pm_to_nm(self, config):
        """Test conversion from pM to nM (divide by 1000)."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_value": [1000.0],
                "standard_units": ["pM"],
            }
        )

        result = extractor.standardize_units(activities, target_unit="nM")

        # 1000 pM = 1 nM
        assert result["standardized_value"].tolist() == [1.0]

    def test_standardize_units_mixed(self, config):
        """Test conversion with mixed units."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standard_value": [100.0, 1.0, 0.001, 10000.0],
                "standard_units": ["nM", "uM", "mM", "pM"],
            }
        )

        result = extractor.standardize_units(activities, target_unit="nM")

        expected = [100.0, 1000.0, 1000.0, 10.0]
        np.testing.assert_array_almost_equal(result["standardized_value"].values, expected)

    def test_compute_pchembl(self, config):
        """Test pChEMBL value computation."""
        extractor = BioactivityExtractor(config)

        # 100 nM = 1e-7 M, pChEMBL = -log10(1e-7) = 7
        activities = pd.DataFrame(
            {
                "standardized_value": [100.0, 1.0, 10000.0],
            }
        )

        result = extractor.compute_pchembl(activities, value_col="standardized_value")

        # pChEMBL = -log10(nM) + 9 = -log10(100) + 9 = -2 + 9 = 7
        expected = [7.0, 9.0, 5.0]
        np.testing.assert_array_almost_equal(result["pchembl_computed"].values, expected)

    def test_compute_pchembl_handles_zero(self, config):
        """Test that zero activity values result in NaN pChEMBL."""
        extractor = BioactivityExtractor(config)

        activities = pd.DataFrame(
            {
                "standardized_value": [0.0, 100.0],
            }
        )

        result = extractor.compute_pchembl(activities, value_col="standardized_value")

        assert np.isnan(result.iloc[0]["pchembl_computed"])
        assert not np.isnan(result.iloc[1]["pchembl_computed"])

    def test_aggregate_by_compound_target(self, config):
        """Test aggregation of multiple measurements."""
        extractor = BioactivityExtractor(config)

        # Multiple measurements for same compound-target pair
        activities = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL1", "CHEMBL2"],
                "uniprot_id": ["P12345", "P12345", "P12345", "P12345"],
                "standard_type": ["IC50", "IC50", "IC50", "IC50"],
                "standardized_value": [100.0, 200.0, 150.0, 500.0],
            }
        )

        result = extractor.aggregate_by_compound_target(activities)

        # Should have 2 rows (2 unique compound-target pairs)
        assert len(result) == 2

        # Check aggregation columns exist
        assert "standardized_value_median" in result.columns
        assert "standardized_value_mean" in result.columns
        assert "standardized_value_count" in result.columns

    def test_extract_for_linked_data_column_renaming(self, config, sample_linked_data):
        """Test that columns are properly renamed in final output."""
        extractor = BioactivityExtractor(config)

        result = extractor.extract_for_linked_data(sample_linked_data)

        # Check key column renames
        if "molecule_chembl_id" in sample_linked_data.columns:
            assert "chembl_id" in result.columns

        if "canonical_smiles" in sample_linked_data.columns:
            assert "smiles" in result.columns

        if "standard_inchi_key" in sample_linked_data.columns:
            assert "inchikey" in result.columns

    def test_extract_for_linked_data_removes_duplicates(self, config):
        """Test that duplicate rows are removed."""
        extractor = BioactivityExtractor(config)

        # Create data with duplicates
        linked = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"],
                "canonical_smiles": ["CCO", "CCO", "CCN"],
                "standard_inchi_key": ["LFQSC-AAAAAA-N", "LFQSC-AAAAAA-N", "OTHER-BBBBBB-N"],
                "uniprot_id": ["P12345", "P12345", "P12345"],
                "standard_type": ["IC50", "IC50", "IC50"],
                "standard_value": [100, 100, 200],
            }
        )

        result = extractor.extract_for_linked_data(linked)

        # Should have 2 unique records
        assert len(result) == 2

    def test_save_output_parquet(self, config, temp_dir):
        """Test saving output as parquet."""
        config.paths.base_dir = temp_dir
        config.output.format = "parquet"
        extractor = BioactivityExtractor(config)

        df = pd.DataFrame(
            {
                "chembl_id": ["CHEMBL1", "CHEMBL2"],
                "smiles": ["CCO", "CCN"],
            }
        )

        output_path = extractor.save_output(df, name="test_output")

        assert output_path.exists()
        assert output_path.suffix == ".parquet"

        # Verify data can be read back
        loaded = pd.read_parquet(output_path)
        assert len(loaded) == 2

    def test_save_output_csv(self, config, temp_dir):
        """Test saving output as CSV."""
        config.paths.base_dir = temp_dir
        config.output.format = "csv"
        extractor = BioactivityExtractor(config)

        df = pd.DataFrame(
            {
                "chembl_id": ["CHEMBL1", "CHEMBL2"],
                "smiles": ["CCO", "CCN"],
            }
        )

        output_path = extractor.save_output(df, name="test_output")

        assert output_path.exists()
        assert output_path.suffix == ".csv"
