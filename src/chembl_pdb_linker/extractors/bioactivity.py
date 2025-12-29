"""Bioactivity data extraction and processing."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class BioactivityExtractor:
    """Extract and process bioactivity data from ChEMBL."""

    # Standard conversion factors to nM
    UNIT_CONVERSIONS = {
        "nM": 1.0,
        "uM": 1000.0,
        "mM": 1000000.0,
        "pM": 0.001,
        "fM": 0.000001,
        "M": 1000000000.0,
    }

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.intermediate_dir = self.paths["intermediate_dir"]
        self.output_dir = self.paths["output_dir"]

    def filter_activities(
        self,
        activities: pd.DataFrame,
        activity_types: Optional[list[str]] = None,
        units: Optional[list[str]] = None,
        min_confidence: int = 9,
    ) -> pd.DataFrame:
        """Filter activities by type, units, and confidence score.

        Args:
            activities: Raw activity DataFrame
            activity_types: List of activity types to include
            units: List of units to include
            min_confidence: Minimum confidence score

        Returns:
            Filtered DataFrame
        """
        logger.info(f"Filtering activities (starting with {len(activities)} records)")

        df = activities.copy()

        # Filter by activity type
        if activity_types is None:
            activity_types = self.config.chembl.activity_types

        if "standard_type" in df.columns:
            df = df[df["standard_type"].isin(activity_types)]
            logger.info(f"After activity type filter: {len(df)} records")

        # Filter by units
        if units is None:
            units = self.config.chembl.standard_units

        if "standard_units" in df.columns:
            df = df[df["standard_units"].isin(units)]
            logger.info(f"After units filter: {len(df)} records")

        # Filter by confidence score
        if "confidence_score" in df.columns:
            df = df[df["confidence_score"] >= min_confidence]
            logger.info(f"After confidence filter: {len(df)} records")

        # Remove null values
        if "standard_value" in df.columns:
            df = df[df["standard_value"].notna()]
            logger.info(f"After removing null values: {len(df)} records")

        return df

    def standardize_units(
        self,
        activities: pd.DataFrame,
        target_unit: str = "nM",
    ) -> pd.DataFrame:
        """Convert all activity values to a standard unit.

        Args:
            activities: Activities DataFrame with standard_value and standard_units
            target_unit: Target unit for standardization

        Returns:
            DataFrame with standardized values
        """
        logger.info(f"Standardizing activity values to {target_unit}")

        df = activities.copy()

        if "standard_units" not in df.columns or "standard_value" not in df.columns:
            logger.warning("Missing standard_units or standard_value columns")
            return df

        # Create standardized value column
        df["standardized_value"] = df["standard_value"]
        df["standardized_unit"] = target_unit

        # Convert based on original unit
        for unit, factor in self.UNIT_CONVERSIONS.items():
            mask = df["standard_units"] == unit
            target_factor = self.UNIT_CONVERSIONS.get(target_unit, 1.0)
            df.loc[mask, "standardized_value"] = (
                df.loc[mask, "standard_value"] * factor / target_factor
            )

        return df

    def aggregate_by_compound_target(
        self,
        activities: pd.DataFrame,
        value_col: str = "standardized_value",
    ) -> pd.DataFrame:
        """Aggregate multiple measurements for the same compound-target pair.

        Args:
            activities: Activities DataFrame
            value_col: Column containing the activity value

        Returns:
            Aggregated DataFrame with median, mean, std, and count
        """
        logger.info("Aggregating activities by compound-target pair")

        # Define grouping columns
        group_cols = []
        for col in [
            "molecule_chembl_id",
            "chembl_id",
            "target_chembl_id",
            "uniprot_id",
            "standard_type",
        ]:
            if col in activities.columns:
                group_cols.append(col)

        if len(group_cols) < 2:
            logger.warning("Insufficient columns for aggregation")
            return activities

        # Aggregate
        agg_dict = {
            value_col: ["median", "mean", "std", "count"],
        }

        # Add other columns to preserve
        preserve_cols = [
            "canonical_smiles",
            "standard_inchi_key",
            "target_name",
            "pchembl_value",
        ]
        for col in preserve_cols:
            if col in activities.columns:
                agg_dict[col] = "first"

        aggregated = activities.groupby(group_cols).agg(agg_dict).reset_index()

        # Flatten column names
        aggregated.columns = [
            "_".join(col).strip("_") if isinstance(col, tuple) else col
            for col in aggregated.columns
        ]

        logger.info(f"Aggregated to {len(aggregated)} compound-target pairs")

        return aggregated

    def compute_pchembl(
        self,
        activities: pd.DataFrame,
        value_col: str = "standardized_value",
    ) -> pd.DataFrame:
        """Compute pChEMBL values (-log10 of activity in M).

        Args:
            activities: Activities DataFrame with values in nM
            value_col: Column containing activity value in nM

        Returns:
            DataFrame with computed pchembl_computed column
        """
        import numpy as np

        df = activities.copy()

        if value_col not in df.columns:
            logger.warning(f"Column {value_col} not found")
            return df

        # Convert nM to M and compute -log10
        # pChEMBL = -log10(value_in_M) = -log10(value_in_nM * 1e-9)
        # = -log10(value_in_nM) - log10(1e-9) = -log10(value_in_nM) + 9
        df["pchembl_computed"] = -np.log10(df[value_col].replace(0, np.nan)) + 9

        return df

    def extract_for_linked_data(
        self,
        linked_data: pd.DataFrame,
    ) -> pd.DataFrame:
        """Extract and format bioactivity data for the final linked dataset.

        Args:
            linked_data: DataFrame with linked ChEMBL-PDB data

        Returns:
            Cleaned DataFrame ready for output
        """
        logger.info("Extracting bioactivity data for final output")

        df = linked_data.copy()

        # Select and rename relevant columns
        column_mapping = {
            "molecule_chembl_id": "chembl_id",
            "chembl_id": "chembl_id",
            "canonical_smiles": "smiles",
            "standard_inchi_key": "inchikey",
            "uniprot_id": "uniprot_id",
            "target_name": "target_name",
            "target_chembl_id": "target_chembl_id",
            "standard_type": "activity_type",
            "standard_value": "activity_value",
            "standardized_value": "activity_value_nm",
            "standard_units": "activity_unit",
            "pchembl_value": "pchembl",
            "pchembl_computed": "pchembl_computed",
            "assay_chembl_id": "assay_chembl_id",
            "confidence_score": "confidence_score",
            "pdb_id": "pdb_id",
            "pdb_ids": "pdb_ids",
            "ligand_code": "pdb_ligand_code",
        }

        # Rename columns that exist
        rename_dict = {k: v for k, v in column_mapping.items() if k in df.columns}
        df = df.rename(columns=rename_dict)

        # Remove duplicates if any
        df = df.drop_duplicates()

        logger.info(f"Extracted {len(df)} records for final output")

        return df

    def save_output(
        self,
        df: pd.DataFrame,
        name: str = "bioactivity_pdb_linked",
    ) -> Path:
        """Save the final output dataset.

        Args:
            df: DataFrame to save
            name: Name for the output file

        Returns:
            Path to saved file
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)

        output_format = self.config.output.format
        compression = self.config.output.compression

        if output_format == "parquet":
            output_path = self.output_dir / f"{name}.parquet"
            df.to_parquet(output_path, compression=compression, index=False)
        elif output_format == "csv":
            output_path = self.output_dir / f"{name}.csv"
            df.to_csv(output_path, index=False)
        else:
            output_path = self.output_dir / f"{name}.parquet"
            df.to_parquet(output_path, compression=compression, index=False)

        logger.info(f"Saved final dataset to {output_path}")
        return output_path
