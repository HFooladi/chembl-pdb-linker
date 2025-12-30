"""UniProt-based protein-level linking between ChEMBL and PDB."""

import logging
from pathlib import Path

import pandas as pd

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class UniProtLinker:
    """Link ChEMBL targets to PDB structures via UniProt IDs."""

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.intermediate_dir = self.paths["intermediate_dir"]

    def link(
        self,
        chembl_activities: pd.DataFrame,
        pdb_sifts: pd.DataFrame,
    ) -> pd.DataFrame:
        """Link ChEMBL activities to PDB structures via UniProt ID.

        Args:
            chembl_activities: ChEMBL activities with UniProt IDs
            pdb_sifts: SIFTS PDB-UniProt mappings

        Returns:
            DataFrame with linked ChEMBL-PDB data
        """
        logger.info("Linking ChEMBL to PDB via UniProt IDs")

        # Ensure we have UniProt IDs in both dataframes
        if "uniprot_id" not in chembl_activities.columns:
            raise ValueError("chembl_activities must contain 'uniprot_id' column")
        if "uniprot_id" not in pdb_sifts.columns:
            raise ValueError("pdb_sifts must contain 'uniprot_id' column")

        # Get unique UniProt IDs from ChEMBL
        chembl_uniprots = set(chembl_activities["uniprot_id"].dropna().unique())
        logger.info(f"ChEMBL has {len(chembl_uniprots)} unique UniProt IDs")

        # Get unique UniProt IDs from PDB
        pdb_uniprots = set(pdb_sifts["uniprot_id"].dropna().unique())
        logger.info(f"PDB has {len(pdb_uniprots)} unique UniProt IDs")

        # Find intersection
        common_uniprots = chembl_uniprots & pdb_uniprots
        logger.info(f"Found {len(common_uniprots)} UniProt IDs in common")

        # Filter both datasets to common UniProt IDs
        chembl_filtered = chembl_activities[
            chembl_activities["uniprot_id"].isin(list(common_uniprots))
        ].copy()
        pdb_filtered = pdb_sifts[pdb_sifts["uniprot_id"].isin(list(common_uniprots))].copy()

        # Aggregate PDB info per UniProt ID (one protein can have multiple structures)
        pdb_agg = (
            pdb_filtered.groupby("uniprot_id")
            .agg(
                {
                    "pdb_id": lambda x: list(set(x)),
                    "chain_id": lambda x: list(set(x)),
                }
            )
            .reset_index()
        )
        pdb_agg = pdb_agg.rename(columns={"pdb_id": "pdb_ids", "chain_id": "chain_ids"})
        pdb_agg["pdb_count"] = pdb_agg["pdb_ids"].apply(len)

        # Merge ChEMBL activities with PDB info
        linked = chembl_filtered.merge(pdb_agg, on="uniprot_id", how="inner")

        logger.info(
            f"Linked {len(linked)} ChEMBL activities to PDB "
            f"(from {len(chembl_activities)} original activities)"
        )

        return linked

    def get_linkage_stats(
        self,
        chembl_activities: pd.DataFrame,
        pdb_sifts: pd.DataFrame,
    ) -> dict:
        """Get statistics about the linkage between ChEMBL and PDB.

        Args:
            chembl_activities: ChEMBL activities with UniProt IDs
            pdb_sifts: SIFTS PDB-UniProt mappings

        Returns:
            Dictionary with linkage statistics
        """
        chembl_uniprots = set(chembl_activities["uniprot_id"].dropna().unique())
        pdb_uniprots = set(pdb_sifts["uniprot_id"].dropna().unique())
        common_uniprots = chembl_uniprots & pdb_uniprots

        chembl_only = chembl_uniprots - pdb_uniprots
        pdb_only = pdb_uniprots - chembl_uniprots

        # Count activities that can be linked
        linkable_activities = chembl_activities[
            chembl_activities["uniprot_id"].isin(list(common_uniprots))
        ]

        return {
            "chembl_uniprot_count": len(chembl_uniprots),
            "pdb_uniprot_count": len(pdb_uniprots),
            "common_uniprot_count": len(common_uniprots),
            "chembl_only_uniprot_count": len(chembl_only),
            "pdb_only_uniprot_count": len(pdb_only),
            "total_activities": len(chembl_activities),
            "linkable_activities": len(linkable_activities),
            "linkage_rate": len(linkable_activities) / len(chembl_activities)
            if len(chembl_activities) > 0
            else 0,
        }

    def save_intermediate(self, df: pd.DataFrame, name: str) -> Path:
        """Save intermediate linked data.

        Args:
            df: DataFrame to save
            name: Name for the output file

        Returns:
            Path to saved file
        """
        self.intermediate_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.intermediate_dir / f"linked_{name}.parquet"
        df.to_parquet(output_path, compression="snappy")
        logger.info(f"Saved {name} to {output_path}")
        return output_path
