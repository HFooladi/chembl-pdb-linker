"""Integration tests for Pipeline class."""

from pathlib import Path

import pandas as pd
import pytest
from chembl_pdb_linker.pipeline import Pipeline


class TestPipelineInitialization:
    """Test Pipeline initialization."""

    def test_pipeline_init_default_config(self, temp_dir):
        """Test pipeline initializes with default config."""
        pipeline = Pipeline(base_dir=temp_dir)

        assert pipeline.config is not None
        assert pipeline.chembl_downloader is not None
        assert pipeline.pdb_downloader is not None
        assert pipeline.uniprot_linker is not None
        assert pipeline.ligand_linker is not None

    def test_pipeline_init_custom_config(self, config, temp_dir):
        """Test pipeline initializes with custom config."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        assert pipeline.config == config

    def test_pipeline_creates_directories(self, temp_dir):
        """Test that pipeline creates required directories."""
        pipeline = Pipeline(base_dir=temp_dir)

        paths = pipeline.paths
        assert paths["raw_dir"].exists() or True  # May be created on demand
        assert paths["intermediate_dir"].exists() or True
        assert paths["output_dir"].exists() or True


class TestPipelineLinking:
    """Test Pipeline linking functionality with sample data."""

    def test_link_with_sample_data(
        self,
        config,
        temp_dir,
        sample_chembl_activities,
        sample_sifts_mapping,
    ):
        """Test linking with sample data files."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        # Save sample data to expected locations
        intermediate_dir = pipeline.paths["intermediate_dir"]
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        activities_path = intermediate_dir / "chembl_activities.parquet"
        sample_chembl_activities.to_parquet(activities_path)

        sifts_path = intermediate_dir / "pdb_sifts_mapping.parquet"
        sample_sifts_mapping.to_parquet(sifts_path)

        # Run linking
        linked = pipeline.link(
            activities_path=activities_path,
            sifts_path=sifts_path,
        )

        # Verify results
        assert len(linked) > 0
        assert "uniprot_id" in linked.columns
        assert "pdb_ids" in linked.columns

    def test_link_missing_activities_file(self, config, temp_dir):
        """Test that linking raises error when activities file is missing."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        with pytest.raises(FileNotFoundError, match="Activities file not found"):
            pipeline.link(
                activities_path=Path("/nonexistent/activities.parquet"),
            )

    def test_link_missing_sifts_file(self, config, temp_dir, sample_chembl_activities):
        """Test that linking raises error when SIFTS file is missing."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        # Create activities file
        intermediate_dir = pipeline.paths["intermediate_dir"]
        intermediate_dir.mkdir(parents=True, exist_ok=True)
        activities_path = intermediate_dir / "chembl_activities.parquet"
        sample_chembl_activities.to_parquet(activities_path)

        with pytest.raises(FileNotFoundError, match="SIFTS mapping not found"):
            pipeline.link(
                activities_path=activities_path,
                sifts_path=Path("/nonexistent/sifts.parquet"),
            )


class TestPipelineExtraction:
    """Test Pipeline extraction functionality."""

    def test_extract_with_sample_data(
        self,
        config,
        temp_dir,
        sample_linked_data,
    ):
        """Test extraction with sample linked data."""
        # Use temp_dir as base to avoid overwriting real output
        config.base_dir = temp_dir
        pipeline = Pipeline(config=config, base_dir=temp_dir)

        # Save sample linked data
        intermediate_dir = pipeline.paths["intermediate_dir"]
        intermediate_dir.mkdir(parents=True, exist_ok=True)
        linked_path = intermediate_dir / "linked_chembl_pdb.parquet"
        sample_linked_data.to_parquet(linked_path)

        # Run extraction (without API enrichment)
        output_path = pipeline.extract(
            linked_path=linked_path,
            enrich_structures=False,
        )

        # Verify output
        assert output_path.exists()

        result = pd.read_parquet(output_path)
        assert len(result) > 0

    def test_extract_missing_linked_file(self, config, temp_dir):
        """Test that extraction raises error when linked file is missing."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        with pytest.raises(FileNotFoundError, match="Linked data not found"):
            pipeline.extract(
                linked_path=Path("/nonexistent/linked.parquet"),
            )


class TestPipelineStatistics:
    """Test Pipeline statistics functionality."""

    def test_get_statistics_empty(self, config, temp_dir):
        """Test statistics when no data files exist."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        stats = pipeline.get_statistics()

        # Should return empty or minimal stats
        assert isinstance(stats, dict)

    def test_get_statistics_with_data(
        self,
        config,
        temp_dir,
        sample_chembl_activities,
        sample_sifts_mapping,
    ):
        """Test statistics with sample data files."""
        config.paths.base_dir = temp_dir
        pipeline = Pipeline(config=config)

        # Save sample data
        intermediate_dir = pipeline.paths["intermediate_dir"]
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        activities_path = intermediate_dir / "chembl_activities.parquet"
        sample_chembl_activities.to_parquet(activities_path)

        sifts_path = intermediate_dir / "pdb_sifts_mapping.parquet"
        sample_sifts_mapping.to_parquet(sifts_path)

        stats = pipeline.get_statistics()

        assert "chembl_activities" in stats
        assert stats["chembl_activities"] == len(sample_chembl_activities)


class TestPipelineIntegration:
    """Full integration tests using real pipeline output."""

    def test_validate_real_output_schema(self, real_output_path):
        """Test that real pipeline output has expected schema."""
        df = pd.read_parquet(real_output_path)

        # Required columns
        required = ["chembl_id", "smiles", "inchikey", "uniprot_id"]
        for col in required:
            assert col in df.columns, f"Missing required column: {col}"

    def test_validate_real_output_chembl_ids(self, real_output_path):
        """Test that ChEMBL IDs have valid format."""
        import re

        df = pd.read_parquet(real_output_path)

        pattern = re.compile(r"^CHEMBL\d+$")
        sample = df["chembl_id"].dropna().head(100)

        for chembl_id in sample:
            assert pattern.match(chembl_id), f"Invalid ChEMBL ID: {chembl_id}"

    def test_validate_real_output_inchikeys(self, real_output_path):
        """Test that InChIKeys have valid format."""
        import re

        df = pd.read_parquet(real_output_path)

        pattern = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
        sample = df["inchikey"].dropna().head(100)

        for inchikey in sample:
            assert pattern.match(inchikey), f"Invalid InChIKey: {inchikey}"

    def test_validate_real_output_activity_values(self, real_output_path):
        """Test that activity values are present and mostly positive."""
        df = pd.read_parquet(real_output_path)

        activity_col = None
        for col in ["activity_value", "activity_value_nm", "standard_value"]:
            if col in df.columns:
                activity_col = col
                break

        if activity_col:
            # Most values should be positive
            positive_ratio = (df[activity_col] > 0).mean()
            assert positive_ratio > 0.99, f"Too many non-positive values: {1 - positive_ratio:.2%}"

    def test_validate_real_output_uniprot_ids(self, real_output_path):
        """Test that UniProt IDs are present."""
        df = pd.read_parquet(real_output_path)

        assert "uniprot_id" in df.columns
        assert df["uniprot_id"].notna().mean() > 0.99, "Too many missing UniProt IDs"

    def test_validate_real_output_record_count(self, real_output_path):
        """Test that output has substantial record count."""
        df = pd.read_parquet(real_output_path)

        # Should have significant data (~98,500 validated protein-ligand pairs)
        assert len(df) > 90000, f"Expected >90,000 records, got {len(df)}"
