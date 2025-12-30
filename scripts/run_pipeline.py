#!/usr/bin/env python3
"""Script to run the ChEMBL-PDB linking pipeline."""

import logging
import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chembl_pdb_linker.config import Config  # noqa: E402
from chembl_pdb_linker.pipeline import Pipeline  # noqa: E402

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def main():
    """Run the complete pipeline."""
    # Get base directory (project root)
    base_dir = Path(__file__).parent.parent

    # Load config
    config = Config.default(base_dir)

    # Create and run pipeline
    pipeline = Pipeline(config)

    print("=" * 60)
    print("ChEMBL-PDB Linker Pipeline")
    print("=" * 60)

    # Run complete pipeline
    output_path = pipeline.run(
        chembl_version="latest",
        enrich_structures=True,
    )

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Output: {output_path}")
    print("=" * 60)

    # Show statistics
    stats = pipeline.get_statistics()
    print("\nDataset Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value:,}" if isinstance(value, int) else f"  {key}: {value}")


if __name__ == "__main__":
    main()
