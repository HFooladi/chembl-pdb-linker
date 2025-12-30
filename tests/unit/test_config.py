"""Tests for configuration module."""

from pathlib import Path

import yaml
from chembl_pdb_linker.config import Config


class TestConfig:
    """Test cases for Config class."""

    def test_default_config_loads(self):
        """Test that default configuration loads without errors."""
        config = Config()

        # Check default values exist
        assert config.chembl is not None
        assert config.pdb is not None
        assert config.linking is not None
        assert config.output is not None

    def test_chembl_defaults(self):
        """Test ChEMBL configuration defaults."""
        config = Config()

        # ChEMBL config has activity_types, confidence_score, etc.
        assert "IC50" in config.chembl.activity_types
        assert "Ki" in config.chembl.activity_types
        assert config.chembl.confidence_score >= 0
        assert "nM" in config.chembl.standard_units
        assert config.chembl.ftp_base is not None

    def test_pdb_defaults(self):
        """Test PDB configuration defaults."""
        config = Config()

        assert config.pdb.max_resolution > 0
        assert config.pdb.api_base is not None
        assert "http" in config.pdb.api_base.lower()

    def test_linking_defaults(self):
        """Test linking configuration defaults."""
        config = Config()

        assert isinstance(config.linking.inchikey_connectivity_only, bool)

    def test_output_defaults(self):
        """Test output configuration defaults."""
        config = Config()

        assert config.output.format in ["parquet", "csv"]
        assert config.output.compression is not None

    def test_get_paths(self):
        """Test path generation."""
        config = Config()
        paths = config.get_paths()

        assert "raw_dir" in paths
        assert "intermediate_dir" in paths
        assert "output_dir" in paths

        # All should be Path objects
        for _key, path in paths.items():
            assert isinstance(path, Path)

    def test_get_paths_with_custom_base_dir(self):
        """Test path generation with custom base directory."""
        config = Config(base_dir=Path("/custom/path"))
        paths = config.get_paths()

        assert "/custom/path" in str(paths["raw_dir"])
        assert "/custom/path" in str(paths["intermediate_dir"])
        assert "/custom/path" in str(paths["output_dir"])

    def test_load_from_yaml(self, temp_dir):
        """Test loading configuration from YAML file."""
        # Create a test YAML config
        yaml_content = {
            "chembl": {
                "activity_types": ["IC50", "Ki"],
                "confidence_score": 8,
            },
            "pdb": {
                "max_resolution": 2.0,
            },
        }

        config_path = temp_dir / "test_config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(yaml_content, f)

        config = Config.from_yaml(config_path, base_dir=temp_dir)

        assert config.chembl.confidence_score == 8
        assert config.pdb.max_resolution == 2.0

    def test_api_config_defaults(self):
        """Test API configuration defaults."""
        config = Config()

        assert config.api.rate_limit > 0
        assert config.api.timeout > 0

    def test_paths_are_relative_to_base(self):
        """Test that data paths are relative to base directory."""
        config = Config(base_dir=Path("/test/base"))
        paths = config.get_paths()

        assert str(paths["raw_dir"]).startswith("/test/base")

    def test_config_immutable_defaults(self):
        """Test that modifying one config doesn't affect others."""
        config1 = Config()
        config2 = Config()

        config1.chembl.confidence_score = 5

        # config2 should still have default
        assert config2.chembl.confidence_score == 9
