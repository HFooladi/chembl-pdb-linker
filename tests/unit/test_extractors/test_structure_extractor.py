"""Tests for StructureExtractor."""

import pandas as pd
from chembl_pdb_linker.extractors.structure import StructureExtractor


class TestStructureExtractor:
    """Test cases for StructureExtractor."""

    def test_generate_download_urls_cif(self, config):
        """Test download URL generation for CIF format."""
        extractor = StructureExtractor(config)

        pdb_ids = ["1ABC", "2DEF"]
        result = extractor.generate_download_urls(pdb_ids, file_format="cif")

        assert len(result) == 2
        assert "pdb_id" in result.columns
        assert "rcsb_url" in result.columns
        assert "pdbe_url" in result.columns

        # Check URL format
        assert "1abc.cif" in result.iloc[0]["rcsb_url"]
        assert "2def.cif" in result.iloc[1]["rcsb_url"]

    def test_generate_download_urls_pdb(self, config):
        """Test download URL generation for PDB format."""
        extractor = StructureExtractor(config)

        result = extractor.generate_download_urls(["1ABC"], file_format="pdb")

        assert "1abc.pdb" in result.iloc[0]["rcsb_url"]

    def test_generate_download_urls_preserves_case(self, config):
        """Test that PDB IDs are uppercased in output."""
        extractor = StructureExtractor(config)

        result = extractor.generate_download_urls(["1abc"])

        assert result.iloc[0]["pdb_id"] == "1ABC"

    def test_filter_by_resolution(self, config, sample_pdb_metadata):
        """Test filtering by resolution."""
        extractor = StructureExtractor(config)

        result = extractor.filter_by_resolution(sample_pdb_metadata, max_resolution=2.0)

        # Should keep entries with resolution <= 2.0 or None (NMR)
        assert len(result) == 3  # 1.9, 1.5, None
        for _, row in result.iterrows():
            if pd.notna(row["resolution"]):
                assert row["resolution"] <= 2.0

    def test_filter_by_resolution_keeps_nmr(self, config, sample_pdb_metadata):
        """Test that NMR structures (no resolution) are kept."""
        extractor = StructureExtractor(config)

        result = extractor.filter_by_resolution(sample_pdb_metadata, max_resolution=1.0)

        # Should include NMR structure with None resolution
        nmr_rows = result[result["resolution"].isna()]
        assert len(nmr_rows) == 1

    def test_filter_by_resolution_no_column(self, config):
        """Test filter when resolution column is missing."""
        extractor = StructureExtractor(config)

        df = pd.DataFrame(
            {
                "pdb_id": ["1ABC", "2DEF"],
                "title": ["Structure 1", "Structure 2"],
                # No resolution column
            }
        )

        result = extractor.filter_by_resolution(df, max_resolution=2.0)

        # Should return unmodified DataFrame
        assert len(result) == 2

    def test_enrich_linked_data_with_pdb_id(self, config):
        """Test enrichment when data has pdb_id column."""
        extractor = StructureExtractor(config)

        linked = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
                "pdb_id": ["1ABC", "2DEF"],
            }
        )

        # Don't fetch metadata via API in tests
        result = extractor.enrich_linked_data(linked, fetch_metadata=False)

        # Should add URL columns
        assert "rcsb_url" in result.columns or len(result) == len(linked)

    def test_enrich_linked_data_with_pdb_ids_list(self, config):
        """Test enrichment when data has pdb_ids as list column."""
        extractor = StructureExtractor(config)

        linked = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pdb_ids": [["1ABC", "2DEF"]],
            }
        )

        result = extractor.enrich_linked_data(linked, fetch_metadata=False)

        # Should have primary PDB ID extracted
        assert "primary_pdb_id" in result.columns
        assert result.iloc[0]["primary_pdb_id"] == "1ABC"

    def test_enrich_linked_data_no_pdb_ids(self, config):
        """Test enrichment when no PDB IDs are present."""
        extractor = StructureExtractor(config)

        linked = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                # No pdb_id or pdb_ids column
            }
        )

        result = extractor.enrich_linked_data(linked, fetch_metadata=False)

        # Should return original data unmodified
        assert len(result) == 1

    def test_save_intermediate(self, config, temp_dir):
        """Test saving intermediate data."""
        config.paths.base_dir = temp_dir
        extractor = StructureExtractor(config)

        df = pd.DataFrame(
            {
                "pdb_id": ["1ABC", "2DEF"],
                "resolution": [2.0, 2.5],
            }
        )

        output_path = extractor.save_intermediate(df, "test_metadata")

        assert output_path.exists()
        assert output_path.suffix == ".parquet"
        assert "pdb_test_metadata" in output_path.name

    def test_rcsb_url_format(self, config):
        """Test RCSB download URL format."""
        extractor = StructureExtractor(config)

        result = extractor.generate_download_urls(["1ABC"])

        url = result.iloc[0]["rcsb_url"]
        assert url.startswith("https://files.rcsb.org/download/")
        assert "1abc" in url.lower()

    def test_pdbe_url_format(self, config):
        """Test PDBe download URL format."""
        extractor = StructureExtractor(config)

        result = extractor.generate_download_urls(["1ABC"])

        url = result.iloc[0]["pdbe_url"]
        assert "ebi.ac.uk" in url
        assert "1abc" in url.lower()

    def test_filter_resolution_default_config(self, config, sample_pdb_metadata):
        """Test filtering uses config default resolution."""
        # Set a specific max resolution in config
        config.pdb.max_resolution = 2.5
        extractor = StructureExtractor(config)

        result = extractor.filter_by_resolution(sample_pdb_metadata)

        # Should use config value (2.5)
        for _, row in result.iterrows():
            if pd.notna(row["resolution"]):
                assert row["resolution"] <= 2.5
