#!/usr/bin/env python3
"""Validate the pipeline output for data quality and schema correctness."""

import re
import sys
from pathlib import Path

import pandas as pd

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class OutputValidator:
    """Validate pipeline output data."""

    # Regex patterns for validation
    CHEMBL_ID_PATTERN = re.compile(r"^CHEMBL\d+$")
    INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
    UNIPROT_PATTERN = re.compile(
        r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4}[0-9]$"
    )
    PDB_ID_PATTERN = re.compile(r"^[0-9A-Za-z]{4}$")
    ACTIVITY_TYPES = {"IC50", "Ki", "Kd", "EC50"}

    def __init__(self, output_path: Path):
        """Initialize validator with output file path."""
        self.output_path = output_path
        self.df = pd.read_parquet(output_path)
        self.errors = []
        self.warnings = []

    def validate_schema(self) -> bool:
        """Check required columns exist."""
        required_cols = [
            "chembl_id",
            "smiles",
            "inchikey",
            "uniprot_id",
            "activity_type",
            "activity_value",
        ]
        missing = [c for c in required_cols if c not in self.df.columns]
        if missing:
            self.errors.append(f"Missing required columns: {missing}")
            return False
        return True

    def validate_chembl_ids(self) -> bool:
        """Validate ChEMBL ID format."""
        invalid_count = (~self.df["chembl_id"].str.match(self.CHEMBL_ID_PATTERN, na=False)).sum()
        if invalid_count > 0:
            self.errors.append(f"Invalid ChEMBL IDs: {invalid_count} records")
            return False
        return True

    def validate_inchikeys(self) -> bool:
        """Validate InChIKey format (27 chars, hyphenated)."""
        invalid_count = (~self.df["inchikey"].str.match(self.INCHIKEY_PATTERN, na=False)).sum()
        if invalid_count > 0:
            self.warnings.append(f"Non-standard InChIKeys: {invalid_count} records")
        return True

    def validate_activity_types(self) -> bool:
        """Check activity types are in expected set."""
        types = set(self.df["activity_type"].dropna().unique())
        unexpected = types - self.ACTIVITY_TYPES
        if unexpected:
            self.warnings.append(f"Unexpected activity types: {unexpected}")
        return True

    def validate_activity_values(self) -> bool:
        """Validate activity values are positive."""
        negative_count = (self.df["activity_value"] <= 0).sum()
        if negative_count > 0:
            self.warnings.append(f"Non-positive activity values: {negative_count} records")
        return True

    def validate_confidence_scores(self) -> bool:
        """Validate confidence scores are 1-9."""
        if "confidence_score" in self.df.columns:
            invalid_count = (~self.df["confidence_score"].between(1, 9, inclusive="both")).sum()
            if invalid_count > 0:
                self.warnings.append(f"Invalid confidence scores: {invalid_count} records")
        return True

    def compute_statistics(self) -> dict:
        """Compute dataset statistics."""
        stats = {
            "total_records": len(self.df),
            "unique_compounds": self.df["chembl_id"].nunique(),
            "unique_targets": self.df["uniprot_id"].nunique(),
        }

        # Activity type distribution
        if "activity_type" in self.df.columns:
            stats["activity_types"] = self.df["activity_type"].value_counts().to_dict()

        # Activity value statistics
        if "activity_value" in self.df.columns:
            stats["activity_value_stats"] = {
                "min": float(self.df["activity_value"].min()),
                "max": float(self.df["activity_value"].max()),
                "mean": float(self.df["activity_value"].mean()),
                "median": float(self.df["activity_value"].median()),
            }

        # Completeness
        stats["completeness"] = {
            col: float(self.df[col].notna().mean() * 100) for col in self.df.columns
        }

        return stats

    def run_all_validations(self) -> dict:
        """Run all validations and return summary."""
        validations = [
            ("schema", self.validate_schema),
            ("chembl_ids", self.validate_chembl_ids),
            ("inchikeys", self.validate_inchikeys),
            ("activity_types", self.validate_activity_types),
            ("activity_values", self.validate_activity_values),
            ("confidence_scores", self.validate_confidence_scores),
        ]

        results = {}
        for name, func in validations:
            try:
                results[name] = func()
            except Exception as e:
                results[name] = False
                self.errors.append(f"Validation {name} failed: {e}")

        return {
            "passed": len(self.errors) == 0,
            "results": results,
            "errors": self.errors,
            "warnings": self.warnings,
            "record_count": len(self.df),
            "column_count": len(self.df.columns),
        }


def main():
    """Run validation on pipeline output."""
    base_dir = Path(__file__).parent.parent
    output_path = base_dir / "data" / "curated" / "bioactivity_pdb_linked.parquet"

    if not output_path.exists():
        print(f"Output file not found: {output_path}")
        print("Run the pipeline first: chembl-pdb-linker run")
        sys.exit(1)

    print("=" * 60)
    print("OUTPUT VALIDATION")
    print("=" * 60)
    print(f"\nFile: {output_path}")

    validator = OutputValidator(output_path)

    # Run validations
    print("\nRunning validations...")
    results = validator.run_all_validations()

    # Print validation results
    print("\n" + "-" * 60)
    print("VALIDATION RESULTS")
    print("-" * 60)
    for name, passed in results["results"].items():
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")

    # Print errors
    if results["errors"]:
        print("\n" + "-" * 60)
        print("ERRORS")
        print("-" * 60)
        for error in results["errors"]:
            print(f"  - {error}")

    # Print warnings
    if results["warnings"]:
        print("\n" + "-" * 60)
        print("WARNINGS")
        print("-" * 60)
        for warning in results["warnings"]:
            print(f"  - {warning}")

    # Print statistics
    print("\n" + "-" * 60)
    print("STATISTICS")
    print("-" * 60)
    stats = validator.compute_statistics()
    print(f"  Total records: {stats['total_records']:,}")
    print(f"  Unique compounds: {stats['unique_compounds']:,}")
    print(f"  Unique targets: {stats['unique_targets']:,}")

    if "activity_types" in stats:
        print("\n  Activity type distribution:")
        for atype, count in stats["activity_types"].items():
            pct = count / stats["total_records"] * 100
            print(f"    {atype}: {count:,} ({pct:.1f}%)")

    if "activity_value_stats" in stats:
        print("\n  Activity value statistics:")
        avs = stats["activity_value_stats"]
        print(f"    Min: {avs['min']:.2f}")
        print(f"    Max: {avs['max']:.2e}")
        print(f"    Mean: {avs['mean']:.2f}")
        print(f"    Median: {avs['median']:.2f}")

    # Overall result
    print("\n" + "=" * 60)
    if results["passed"]:
        print("VALIDATION PASSED")
    else:
        print("VALIDATION FAILED")
    print("=" * 60)

    return 0 if results["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
