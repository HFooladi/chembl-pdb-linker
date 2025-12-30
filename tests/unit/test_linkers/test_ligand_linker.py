"""Tests for LigandLinker."""

import pandas as pd
import pytest
from chembl_pdb_linker.linkers.ligand import LigandLinker


class TestLigandLinker:
    """Test cases for LigandLinker."""

    def test_link_full_inchikey(self, config, sample_chembl_activities, sample_ligand_data):
        """Test linking with full InChIKey matching."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        result = linker.link(sample_chembl_activities, sample_ligand_data)

        # Should find matches for aspirin and ibuprofen (2 compounds)
        assert len(result) > 0
        assert "ligand_code" in result.columns

    def test_link_connectivity_only(self, config):
        """Test linking with connectivity-only (first 14 chars) InChIKey matching."""
        config.linking.inchikey_connectivity_only = True
        linker = LigandLinker(config)

        # Same connectivity layer, different stereochemistry
        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "standard_inchi_key": ["BSYNRYMUTXBXSQ-ABCDEFGHIJ-N"],  # Different suffix
            }
        )

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["ASP"],
                "ligand_inchikey": ["BSYNRYMUTXBXSQ-ZYXWVUTSRQ-M"],  # Same first 14 chars
            }
        )

        result = linker.link(chembl, pdb)

        # Should match because first 14 characters are the same
        assert len(result) == 1

    def test_link_no_matches(self, config):
        """Test linking when there are no matching InChIKeys."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "standard_inchi_key": ["AAAAAAAAAAAAAA-BBBBBBBBBB-C"],
            }
        )

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["LIG"],
                "ligand_inchikey": ["ZZZZZZZZZZZZZ-YYYYYYYYYY-X"],
            }
        )

        result = linker.link(chembl, pdb)

        # Should return empty DataFrame
        assert len(result) == 0

    def test_link_missing_inchikey_column_chembl(self, config, sample_ligand_data):
        """Test that missing InChIKey column in ChEMBL raises ValueError."""
        linker = LigandLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                # No InChIKey column
            }
        )

        with pytest.raises(ValueError, match="must contain an InChIKey column"):
            linker.link(chembl, sample_ligand_data)

    def test_link_missing_inchikey_column_pdb(self, config, sample_chembl_activities):
        """Test that missing InChIKey column in PDB raises ValueError."""
        linker = LigandLinker(config)

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["LIG"],
                # No InChIKey column
            }
        )

        with pytest.raises(ValueError, match="must contain an InChIKey column"):
            linker.link(sample_chembl_activities, pdb)

    def test_get_inchikey_prefix(self, config):
        """Test InChIKey prefix extraction."""
        linker = LigandLinker(config)

        # Normal InChIKey
        assert linker._get_inchikey_prefix("BSYNRYMUTXBXSQ-UHFFFAOYSA-N") == "BSYNRYMUTXBXSQ"

        # Short string
        assert linker._get_inchikey_prefix("SHORT") == "SHORT"

        # Empty/None
        assert linker._get_inchikey_prefix("") == ""
        assert linker._get_inchikey_prefix(None) == ""

    def test_link_with_activities(self, config, sample_chembl_activities, sample_ligand_data):
        """Test linking activities to specific PDB ligand structures."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        result = linker.link_with_activities(sample_chembl_activities, sample_ligand_data)

        # Should have linked some activities
        if len(result) > 0:
            assert "pdb_id" in result.columns
            assert "ligand_code" in result.columns
            assert "molecule_chembl_id" in result.columns

    def test_get_linkage_stats(self, config, sample_chembl_activities, sample_ligand_data):
        """Test linkage statistics computation."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        stats = linker.get_linkage_stats(sample_chembl_activities, sample_ligand_data)

        # Should have expected keys
        assert "chembl_inchikey_count" in stats
        assert "pdb_inchikey_count" in stats
        assert "common_inchikey_count" in stats
        assert "match_type" in stats
        assert "ligand_linkage_rate" in stats
        assert "has_protein_linking" in stats

        # Match type should be full
        assert stats["match_type"] == "full_inchikey"

    def test_get_linkage_stats_connectivity_only(
        self, config, sample_chembl_activities, sample_ligand_data
    ):
        """Test linkage statistics with connectivity-only matching."""
        config.linking.inchikey_connectivity_only = True
        linker = LigandLinker(config)

        stats = linker.get_linkage_stats(sample_chembl_activities, sample_ligand_data)

        assert stats["match_type"] == "connectivity_only"

    def test_alternative_column_names(self, config):
        """Test that alternative column names are recognized."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        # Use 'inchikey' instead of 'standard_inchi_key'
        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
            }
        )

        # Use 'inchi_key' instead of 'ligand_inchikey'
        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["ASP"],
                "inchi_key": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
            }
        )

        result = linker.link(chembl, pdb)

        # Should work with alternative column names
        assert len(result) == 1

    def test_handles_nan_inchikeys(self, config):
        """Test that NaN InChIKeys are handled gracefully."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        chembl = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "standard_inchi_key": [
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    None,
                    "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                ],
            }
        )

        pdb = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["ASP"],
                "ligand_inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
            }
        )

        result = linker.link(chembl, pdb)

        # Should match CHEMBL1 only (NaN skipped, CHEMBL3 doesn't match)
        assert len(result) == 1

    def test_save_intermediate(
        self, config, sample_chembl_activities, sample_ligand_data, temp_dir
    ):
        """Test saving intermediate linked data."""
        config.paths.base_dir = temp_dir
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        result = linker.link(sample_chembl_activities, sample_ligand_data)
        if len(result) > 0:
            output_path = linker.save_intermediate(result, "test_ligand_linked")

            assert output_path.exists()
            assert output_path.suffix == ".parquet"


class TestProteinLigandPairValidation:
    """Test cases for protein-ligand pair validation in LigandLinker."""

    def test_validate_protein_ligand_pair_valid(self, config):
        """Test validation returns True when pdb_id is in pdb_ids list."""
        linker = LigandLinker(config)

        # PDB ID is in the list
        assert linker._validate_protein_ligand_pair("1ABC", ["1ABC", "2DEF", "3GHI"]) is True
        # Case insensitive
        assert linker._validate_protein_ligand_pair("1abc", ["1ABC", "2DEF"]) is True
        assert linker._validate_protein_ligand_pair("1ABC", ["1abc", "2def"]) is True

    def test_validate_protein_ligand_pair_invalid(self, config):
        """Test validation returns False when pdb_id is not in pdb_ids list."""
        linker = LigandLinker(config)

        # PDB ID is not in the list
        assert linker._validate_protein_ligand_pair("1XYZ", ["1ABC", "2DEF", "3GHI"]) is False
        # Empty list
        assert linker._validate_protein_ligand_pair("1ABC", []) is False
        # None values
        assert linker._validate_protein_ligand_pair(None, ["1ABC"]) is False
        assert linker._validate_protein_ligand_pair("1ABC", None) is False
        # String instead of list
        assert linker._validate_protein_ligand_pair("1ABC", "1ABC") is False

    def test_link_with_activities_validates_protein_ligand_pairs(self, config):
        """Test that link_with_activities validates protein-ligand pairs.

        This is the critical test: when activities have pdb_ids from protein linking,
        only keep ligand matches where the PDB structure contains both protein and ligand.
        """
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        # Activities with protein linking info (pdb_ids column)
        # Protein P23219 has structures in 1PRH and 4U0J
        activities = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL25"],
                "standard_inchi_key": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],  # Aspirin
                "uniprot_id": ["P23219"],
                "pdb_ids": [["1PRH", "4U0J"]],  # PDBs from protein linking
            }
        )

        # PDB ligands - aspirin is in 1PRH (matches protein) and 9ZZZ (doesn't match protein)
        pdb_ligands = pd.DataFrame(
            {
                "pdb_id": ["1PRH", "9ZZZ"],
                "ligand_code": ["ASP", "ASP"],
                "ligand_inchikey": [
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",  # Aspirin in 1PRH (valid)
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",  # Aspirin in 9ZZZ (invalid - no protein)
                ],
            }
        )

        result = linker.link_with_activities(activities, pdb_ligands)

        # Should only return 1PRH match, not 9ZZZ
        assert len(result) == 1
        assert result.iloc[0]["pdb_id"] == "1PRH"

    def test_link_with_activities_no_valid_pairs(self, config):
        """Test when ligand and protein are in different PDB structures."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        # Protein is in 1ABC, 2DEF
        activities = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL25"],
                "standard_inchi_key": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
                "uniprot_id": ["P23219"],
                "pdb_ids": [["1ABC", "2DEF"]],
            }
        )

        # Ligand is in 3GHI, 4JKL (different structures, no overlap)
        pdb_ligands = pd.DataFrame(
            {
                "pdb_id": ["3GHI", "4JKL"],
                "ligand_code": ["ASP", "ASP"],
                "ligand_inchikey": [
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                ],
            }
        )

        result = linker.link_with_activities(activities, pdb_ligands)

        # Should return empty - no PDB structure has both protein and ligand
        assert len(result) == 0

    def test_link_with_activities_without_protein_linking(self, config):
        """Test behavior when no protein linking has been done (no pdb_ids column)."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        # Activities without pdb_ids column (no protein linking done)
        activities = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL25"],
                "standard_inchi_key": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
                "uniprot_id": ["P23219"],
            }
        )

        pdb_ligands = pd.DataFrame(
            {
                "pdb_id": ["1ABC"],
                "ligand_code": ["ASP"],
                "ligand_inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
            }
        )

        result = linker.link_with_activities(activities, pdb_ligands)

        # Should still return matches (but warns about missing protein linking)
        assert len(result) == 1

    def test_link_with_activities_multiple_valid_pairs(self, config):
        """Test with multiple valid protein-ligand pairs."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        # Two activities for same ligand with different proteins/PDBs
        activities = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL25", "CHEMBL25"],
                "standard_inchi_key": [
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                ],
                "uniprot_id": ["P23219", "P35354"],
                "pdb_ids": [["1PRH", "4U0J"], ["3PGH"]],
            }
        )

        # Ligand is in 1PRH and 3PGH
        pdb_ligands = pd.DataFrame(
            {
                "pdb_id": ["1PRH", "3PGH", "9ZZZ"],
                "ligand_code": ["ASP", "ASP", "ASP"],
                "ligand_inchikey": [
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                ],
            }
        )

        result = linker.link_with_activities(activities, pdb_ligands)

        # Should return 2 matches: 1PRH for first activity, 3PGH for second
        assert len(result) == 2
        pdb_ids = set(result["pdb_id"].tolist())
        assert pdb_ids == {"1PRH", "3PGH"}
        # 9ZZZ should not be in results
        assert "9ZZZ" not in pdb_ids

    def test_get_linkage_stats_with_protein_linking(
        self, config, sample_linked_data, sample_ligand_data
    ):
        """Test linkage statistics include validated pair counts."""
        config.linking.inchikey_connectivity_only = False
        linker = LigandLinker(config)

        stats = linker.get_linkage_stats(sample_linked_data, sample_ligand_data)

        # Should have protein linking stats
        assert "has_protein_linking" in stats
        assert stats["has_protein_linking"] is True

        # Should have validated pair stats if there were matches
        if stats.get("common_inchikey_count", 0) > 0:
            assert "activities_with_ligand_match" in stats
            assert "activities_with_validated_pairs" in stats
            assert "validated_pair_rate" in stats
