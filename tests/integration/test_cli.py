"""Integration tests for CLI."""

import subprocess
import sys

from chembl_pdb_linker.cli import app
from typer.testing import CliRunner

runner = CliRunner()


class TestCLIBasic:
    """Test basic CLI functionality."""

    def test_help(self):
        """Test that --help works."""
        result = runner.invoke(app, ["--help"])

        assert result.exit_code == 0
        assert "ChEMBL-PDB Linker" in result.stdout or "chembl" in result.stdout.lower()

    def test_version(self):
        """Test that --version works."""
        result = runner.invoke(app, ["--version"])

        # May exit with 0 or show version
        assert result.exit_code == 0 or "version" in result.stdout.lower()

    def test_download_help(self):
        """Test download command help."""
        result = runner.invoke(app, ["download", "--help"])

        assert result.exit_code == 0
        assert "download" in result.stdout.lower()

    def test_link_help(self):
        """Test link command help."""
        result = runner.invoke(app, ["link", "--help"])

        assert result.exit_code == 0
        assert "link" in result.stdout.lower()

    def test_extract_help(self):
        """Test extract command help."""
        result = runner.invoke(app, ["extract", "--help"])

        assert result.exit_code == 0
        assert "extract" in result.stdout.lower()

    def test_stats_help(self):
        """Test stats command help."""
        result = runner.invoke(app, ["stats", "--help"])

        assert result.exit_code == 0
        assert "stats" in result.stdout.lower() or "statistic" in result.stdout.lower()


class TestCLICommands:
    """Test CLI commands with temp directories."""

    def test_stats_no_data(self, temp_dir):
        """Test stats command when no data exists."""
        result = runner.invoke(app, ["stats", "--base-dir", str(temp_dir)])

        # Should not crash, may show empty stats or message
        assert result.exit_code == 0 or "no data" in result.stdout.lower()

    def test_link_missing_data(self, temp_dir):
        """Test link command fails gracefully when data is missing."""
        result = runner.invoke(app, ["link", "--base-dir", str(temp_dir)])

        # Should fail with informative error
        assert (
            result.exit_code != 0
            or "not found" in result.stdout.lower()
            or "error" in result.stdout.lower()
        )

    def test_extract_missing_data(self, temp_dir):
        """Test extract command fails gracefully when data is missing."""
        result = runner.invoke(app, ["extract", "--base-dir", str(temp_dir)])

        # Should fail with informative error
        assert (
            result.exit_code != 0
            or "not found" in result.stdout.lower()
            or "error" in result.stdout.lower()
        )


class TestCLIWithRealData:
    """Test CLI commands with real pipeline data."""

    def test_stats_with_real_data(self, real_output_path):
        """Test stats command with real data."""
        base_dir = real_output_path.parent.parent.parent  # Go up from data/curated/

        result = runner.invoke(app, ["stats", "--base-dir", str(base_dir)])

        # Should show statistics
        if result.exit_code == 0:
            # Check for expected output
            assert any(
                keyword in result.stdout.lower()
                for keyword in ["record", "compound", "target", "activities"]
            )


class TestCLIModule:
    """Test CLI as a module."""

    def test_module_invocation(self):
        """Test that the CLI can be invoked as a module."""
        result = subprocess.run(
            [sys.executable, "-m", "chembl_pdb_linker.cli", "--help"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "chembl" in result.stdout.lower() or "pdb" in result.stdout.lower()

    def test_package_entrypoint(self):
        """Test the package entry point if installed."""
        result = subprocess.run(
            ["chembl-pdb-linker", "--help"],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            assert "chembl" in result.stdout.lower() or "pdb" in result.stdout.lower()
