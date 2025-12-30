"""Tests for UniProtLinker."""

import pandas as pd
import pytest
from chembl_pdb_linker.linkers.uniprot import UniProtLinker


class TestUniProtLinker:
    """Test cases for UniProtLinker."""

    def test_link_basic(self, config, sample_chembl_activities, sample_sifts_mapping):
        """Test basic linking of ChEMBL activities to PDB structures."""
        linker = UniProtLinker(config)
        result = linker.link(sample_chembl_activities, sample_sifts_mapping)

        # Should have linked records
        assert len(result) > 0

        # Should have required columns
        assert "uniprot_id" in result.columns
        assert "pdb_ids" in result.columns
        assert "pdb_count" in result.columns

        # All records should have PDB IDs
        assert result["pdb_ids"].notna().all()

    def test_link_preserves_activity_columns(
        self, config, sample_chembl_activities, sample_sifts_mapping
    ):
        """Test that original activity columns are preserved after linking."""
        linker = UniProtLinker(config)
        result = linker.link(sample_chembl_activities, sample_sifts_mapping)

        # Original columns should be preserved
        for col in ["molecule_chembl_id", "canonical_smiles", "standard_value", "standard_type"]:
            assert col in result.columns

    def test_link_no_overlap(self, config):
        """Test linking when there are no common UniProt IDs."""
        linker = UniProtLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
                "uniprot_id": ["P00001", "P00002"],
                "standard_value": [100, 200],
            }
        )

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC", "2DEF"],
                "chain_id": ["A", "B"],
                "uniprot_id": ["P99999", "P99998"],  # No overlap
            }
        )

        result = linker.link(chembl, pdb)

        # Should return empty DataFrame
        assert len(result) == 0

    def test_link_missing_uniprot_column_chembl(self, config, sample_sifts_mapping):
        """Test that missing uniprot_id column in ChEMBL raises ValueError."""
        linker = UniProtLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "standard_value": [100],
                # Missing uniprot_id
            }
        )

        with pytest.raises(ValueError, match="chembl_activities must contain 'uniprot_id' column"):
            linker.link(chembl, sample_sifts_mapping)

    def test_link_missing_uniprot_column_pdb(self, config, sample_chembl_activities):
        """Test that missing uniprot_id column in PDB raises ValueError."""
        linker = UniProtLinker(config)

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "chain_id": ["A"],
                # Missing uniprot_id
            }
        )

        with pytest.raises(ValueError, match="pdb_sifts must contain 'uniprot_id' column"):
            linker.link(sample_chembl_activities, pdb)

    def test_link_aggregates_pdb_ids(self, config):
        """Test that multiple PDB structures for same UniProt are aggregated."""
        linker = UniProtLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "uniprot_id": ["P12345"],
                "standard_value": [100],
            }
        )

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC", "2DEF", "3GHI"],
                "chain_id": ["A", "A", "B"],
                "uniprot_id": ["P12345", "P12345", "P12345"],  # Same UniProt, multiple PDBs
            }
        )

        result = linker.link(chembl, pdb)

        assert len(result) == 1
        assert len(result.iloc[0]["pdb_ids"]) == 3
        assert result.iloc[0]["pdb_count"] == 3

    def test_get_linkage_stats(self, config, sample_chembl_activities, sample_sifts_mapping):
        """Test linkage statistics computation."""
        linker = UniProtLinker(config)
        stats = linker.get_linkage_stats(sample_chembl_activities, sample_sifts_mapping)

        # Should have expected keys
        assert "chembl_uniprot_count" in stats
        assert "pdb_uniprot_count" in stats
        assert "common_uniprot_count" in stats
        assert "linkage_rate" in stats

        # Values should be sensible
        assert stats["chembl_uniprot_count"] > 0
        assert stats["pdb_uniprot_count"] > 0
        assert 0 <= stats["linkage_rate"] <= 1

    def test_get_linkage_stats_empty_data(self, config):
        """Test linkage statistics with empty data."""
        linker = UniProtLinker(config)

        chembl = pd.DataFrame({"molecule_chembl_id": [], "uniprot_id": []})
        pdb = pd.DataFrame({"pdb_id": [], "chain_id": [], "uniprot_id": []})

        stats = linker.get_linkage_stats(chembl, pdb)

        assert stats["chembl_uniprot_count"] == 0
        assert stats["linkage_rate"] == 0

    def test_save_intermediate(
        self, config, sample_chembl_activities, sample_sifts_mapping, temp_dir
    ):
        """Test saving intermediate linked data."""
        config.paths.base_dir = temp_dir
        linker = UniProtLinker(config)

        result = linker.link(sample_chembl_activities, sample_sifts_mapping)
        output_path = linker.save_intermediate(result, "test_linked")

        assert output_path.exists()
        assert output_path.suffix == ".parquet"

        # Verify data can be read back
        loaded = pd.read_parquet(output_path)
        assert len(loaded) == len(result)
