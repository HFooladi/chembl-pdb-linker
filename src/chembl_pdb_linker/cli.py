"""Command-line interface for ChEMBL-PDB Linker."""

import logging
import sys
from pathlib import Path
from typing import Optional

import typer

from chembl_pdb_linker import __version__
from chembl_pdb_linker.config import Config
from chembl_pdb_linker.pipeline import Pipeline

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

app = typer.Typer(
    name="chembl-pdb-linker",
    help="Link ChEMBL bioactivity data with PDB structural information.",
    add_completion=False,
)


def get_config(config_path: Optional[Path], base_dir: Optional[Path]) -> Config:
    """Load configuration from file or use defaults."""
    if base_dir is None:
        base_dir = Path.cwd()

    if config_path is not None:
        return Config.from_yaml(config_path, base_dir)
    return Config.default(base_dir)


@app.command()
def download(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to configuration file",
    ),
    base_dir: Optional[Path] = typer.Option(
        None,
        "--base-dir",
        "-d",
        help="Base directory for data files",
    ),
    chembl_version: str = typer.Option(
        "latest",
        "--chembl-version",
        "-v",
        help="ChEMBL version to download (e.g., '34' or 'latest')",
    ),
) -> None:
    """Download ChEMBL and PDB data."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")
    typer.echo("=" * 50)
    typer.echo("Downloading data...")

    cfg = get_config(config, base_dir)
    pipeline = Pipeline(cfg)

    try:
        outputs = pipeline.download(chembl_version=chembl_version)
        typer.echo("\nDownload complete!")
        typer.echo("Downloaded files:")
        for name, path in outputs.items():
            typer.echo(f"  - {name}: {path}")
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from None


@app.command()
def link(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to configuration file",
    ),
    base_dir: Optional[Path] = typer.Option(
        None,
        "--base-dir",
        "-d",
        help="Base directory for data files",
    ),
    activities: Optional[Path] = typer.Option(
        None,
        "--activities",
        "-a",
        help="Path to ChEMBL activities parquet file",
    ),
    sifts: Optional[Path] = typer.Option(
        None,
        "--sifts",
        "-s",
        help="Path to SIFTS mapping parquet file",
    ),
    ligands: Optional[Path] = typer.Option(
        None,
        "--ligands",
        "-l",
        help="Path to PDB ligands parquet file",
    ),
) -> None:
    """Link ChEMBL and PDB data via UniProt and InChIKey."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")
    typer.echo("=" * 50)
    typer.echo("Linking data...")

    cfg = get_config(config, base_dir)
    pipeline = Pipeline(cfg)

    try:
        linked = pipeline.link(
            activities_path=activities,
            sifts_path=sifts,
            ligands_path=ligands,
        )
        typer.echo("\nLinking complete!")
        typer.echo(f"Linked {len(linked)} records")
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        typer.echo("Run 'chembl-pdb-linker download' first.", err=True)
        raise typer.Exit(1) from None
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from None


@app.command()
def extract(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to configuration file",
    ),
    base_dir: Optional[Path] = typer.Option(
        None,
        "--base-dir",
        "-d",
        help="Base directory for data files",
    ),
    linked: Optional[Path] = typer.Option(
        None,
        "--linked",
        "-l",
        help="Path to linked data parquet file",
    ),
    no_enrich: bool = typer.Option(
        False,
        "--no-enrich",
        help="Skip fetching additional structure metadata",
    ),
) -> None:
    """Extract and format the final dataset."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")
    typer.echo("=" * 50)
    typer.echo("Extracting final dataset...")

    cfg = get_config(config, base_dir)
    pipeline = Pipeline(cfg)

    try:
        output_path = pipeline.extract(
            linked_path=linked,
            enrich_structures=not no_enrich,
        )
        typer.echo("\nExtraction complete!")
        typer.echo(f"Output file: {output_path}")
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        typer.echo("Run 'chembl-pdb-linker link' first.", err=True)
        raise typer.Exit(1) from None
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from None


@app.command()
def run(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to configuration file",
    ),
    base_dir: Optional[Path] = typer.Option(
        None,
        "--base-dir",
        "-d",
        help="Base directory for data files",
    ),
    chembl_version: str = typer.Option(
        "latest",
        "--chembl-version",
        "-v",
        help="ChEMBL version to download",
    ),
    no_enrich: bool = typer.Option(
        False,
        "--no-enrich",
        help="Skip fetching additional structure metadata",
    ),
) -> None:
    """Run the complete pipeline (download, link, extract)."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")
    typer.echo("=" * 50)
    typer.echo("Running complete pipeline...")

    cfg = get_config(config, base_dir)
    pipeline = Pipeline(cfg)

    try:
        output_path = pipeline.run(
            chembl_version=chembl_version,
            enrich_structures=not no_enrich,
        )
        typer.echo("\nPipeline complete!")
        typer.echo(f"Output file: {output_path}")
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from None


@app.command()
def stats(
    config: Optional[Path] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to configuration file",
    ),
    base_dir: Optional[Path] = typer.Option(
        None,
        "--base-dir",
        "-d",
        help="Base directory for data files",
    ),
) -> None:
    """Show statistics about the current dataset."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")
    typer.echo("=" * 50)
    typer.echo("Dataset Statistics")
    typer.echo("-" * 50)

    cfg = get_config(config, base_dir)
    pipeline = Pipeline(cfg)

    stats = pipeline.get_statistics()

    if not stats:
        typer.echo("No data found. Run the pipeline first.")
        return

    for key, value in stats.items():
        label = key.replace("_", " ").title()
        if isinstance(value, float):
            typer.echo(f"{label}: {value:.4f}")
        else:
            typer.echo(f"{label}: {value:,}")


@app.command()
def version() -> None:
    """Show version information."""
    typer.echo(f"ChEMBL-PDB Linker v{__version__}")


if __name__ == "__main__":
    app()
