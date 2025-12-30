"""Ligand-based linking between ChEMBL and PDB using InChIKey."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class LigandLinker:
    """Link ChEMBL compounds to PDB ligands via InChIKey matching."""

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.intermediate_dir = self.paths["intermediate_dir"]
        self.connectivity_only = config.linking.inchikey_connectivity_only

    def _get_inchikey_prefix(self, inchikey: str) -> str:
        """Extract the connectivity layer (first 14 chars) from InChIKey.

        Args:
            inchikey: Full InChIKey string

        Returns:
            First 14 characters (connectivity layer)
        """
        if pd.isna(inchikey) or not isinstance(inchikey, str):  # type: ignore[arg-type]
            return ""
        return inchikey[:14] if len(inchikey) >= 14 else inchikey

    def _validate_protein_ligand_pair(
        self, pdb_id: Optional[str], pdb_ids_list: Optional[list]
    ) -> bool:
        """Check if a PDB ID is in the protein's PDB ID list.

        This ensures the PDB structure contains BOTH the target protein AND the ligand.

        Args:
            pdb_id: PDB ID from ligand matching
            pdb_ids_list: List of PDB IDs from protein-level matching

        Returns:
            True if pdb_id is in pdb_ids_list, False otherwise
        """
        if pdb_id is None or pdb_ids_list is None:
            return False
        if not hasattr(pdb_ids_list, "__iter__") or isinstance(pdb_ids_list, str):
            return False
        # Normalize to uppercase for comparison
        pdb_id_upper = pdb_id.upper()
        pdb_ids_upper = [p.upper() if isinstance(p, str) else p for p in pdb_ids_list]
        return pdb_id_upper in pdb_ids_upper

    def link(
        self,
        chembl_compounds: pd.DataFrame,
        pdb_ligands: pd.DataFrame,
    ) -> pd.DataFrame:
        """Link ChEMBL compounds to PDB ligands via InChIKey.

        Args:
            chembl_compounds: ChEMBL compounds with InChIKeys
            pdb_ligands: PDB ligands with InChIKeys

        Returns:
            DataFrame mapping ChEMBL compounds to PDB ligands
        """
        logger.info("Linking ChEMBL compounds to PDB ligands via InChIKey")

        # Validate required columns
        chembl_inchikey_col = None
        for col in ["standard_inchi_key", "inchikey", "standard_inchikey"]:
            if col in chembl_compounds.columns:
                chembl_inchikey_col = col
                break

        if chembl_inchikey_col is None:
            raise ValueError(
                "chembl_compounds must contain an InChIKey column "
                "(standard_inchi_key, inchikey, or standard_inchikey)"
            )

        pdb_inchikey_col = None
        for col in ["ligand_inchikey", "inchikey", "inchi_key"]:
            if col in pdb_ligands.columns:
                pdb_inchikey_col = col
                break

        if pdb_inchikey_col is None:
            raise ValueError(
                "pdb_ligands must contain an InChIKey column "
                "(ligand_inchikey, inchikey, or inchi_key)"
            )

        # Prepare ChEMBL compounds
        chembl_df = chembl_compounds.copy()
        chembl_df["chembl_inchikey"] = chembl_df[chembl_inchikey_col]
        chembl_df = chembl_df.dropna(subset=["chembl_inchikey"])

        # Prepare PDB ligands
        pdb_df = pdb_ligands.copy()
        pdb_df["pdb_inchikey"] = pdb_df[pdb_inchikey_col]
        pdb_df = pdb_df.dropna(subset=["pdb_inchikey"])

        if self.connectivity_only:
            logger.info("Using InChIKey connectivity layer (first 14 chars) for matching")
            chembl_df["match_key"] = chembl_df["chembl_inchikey"].apply(self._get_inchikey_prefix)
            pdb_df["match_key"] = pdb_df["pdb_inchikey"].apply(self._get_inchikey_prefix)
        else:
            logger.info("Using full InChIKey for matching")
            chembl_df["match_key"] = chembl_df["chembl_inchikey"]
            pdb_df["match_key"] = pdb_df["pdb_inchikey"]

        # Find matching InChIKeys
        chembl_keys = set(chembl_df["match_key"].dropna().unique())
        pdb_keys = set(pdb_df["match_key"].dropna().unique())
        common_keys = chembl_keys & pdb_keys

        logger.info(f"ChEMBL has {len(chembl_keys)} unique InChIKeys")
        logger.info(f"PDB has {len(pdb_keys)} unique InChIKeys")
        logger.info(f"Found {len(common_keys)} InChIKeys in common")

        if len(common_keys) == 0:
            logger.warning("No matching InChIKeys found between ChEMBL and PDB")
            return pd.DataFrame()

        # Filter to common keys
        chembl_filtered = chembl_df[chembl_df["match_key"].isin(list(common_keys))]
        pdb_filtered = pdb_df[pdb_df["match_key"].isin(list(common_keys))]

        # Merge on match_key
        linked = chembl_filtered.merge(
            pdb_filtered,
            on="match_key",
            how="inner",
            suffixes=("_chembl", "_pdb"),
        )

        # Clean up columns
        cols_to_keep = [
            col
            for col in linked.columns
            if not col.startswith("match_key") and col != "chembl_inchikey"
        ]
        linked = linked[cols_to_keep]

        logger.info(
            f"Created {len(linked)} ChEMBL-PDB ligand mappings "
            f"(from {len(chembl_compounds)} ChEMBL compounds)"
        )

        return linked  # type: ignore[return-value]

    def link_with_activities(
        self,
        chembl_activities: pd.DataFrame,
        pdb_ligands: pd.DataFrame,
    ) -> pd.DataFrame:
        """Link ChEMBL activities to PDB structures where the ligand matches.

        This combines protein-level and ligand-level matching:
        - Activity must be for a protein with a PDB structure
        - The compound in the activity must match a ligand in that PDB structure
        - The PDB structure must contain BOTH the target protein AND the ligand

        Args:
            chembl_activities: ChEMBL activities with compound InChIKeys and pdb_ids
                              (from protein-level linking)
            pdb_ligands: PDB ligands with InChIKeys and PDB IDs

        Returns:
            DataFrame of activities linked to specific PDB-ligand complexes
            where the structure contains both the protein and ligand
        """
        logger.info("Linking ChEMBL activities to PDB ligand structures")

        # Check if protein-level linking has been done (pdb_ids column present)
        has_protein_linking = "pdb_ids" in chembl_activities.columns

        if not has_protein_linking:
            logger.warning(
                "No 'pdb_ids' column found in activities. "
                "Protein-level linking may not have been performed. "
                "Results may include PDB structures that don't contain the target protein."
            )

        # Identify InChIKey column in activities
        inchikey_col = None
        for col in ["standard_inchi_key", "inchikey", "standard_inchikey"]:
            if col in chembl_activities.columns:
                inchikey_col = col
                break

        if inchikey_col is None:
            raise ValueError("chembl_activities must contain an InChIKey column")

        # Prepare data
        activities = chembl_activities.copy()
        activities["chembl_inchikey"] = activities[inchikey_col]

        pdb_df = pdb_ligands.copy()
        pdb_inchikey_col = None
        for col in ["ligand_inchikey", "inchikey", "inchi_key"]:
            if col in pdb_df.columns:
                pdb_inchikey_col = col
                break

        if pdb_inchikey_col is None:
            raise ValueError("pdb_ligands must contain an InChIKey column")

        pdb_df["pdb_inchikey"] = pdb_df[pdb_inchikey_col]

        # Create match keys
        if self.connectivity_only:
            activities["match_key"] = activities["chembl_inchikey"].apply(self._get_inchikey_prefix)
            pdb_df["match_key"] = pdb_df["pdb_inchikey"].apply(self._get_inchikey_prefix)
        else:
            activities["match_key"] = activities["chembl_inchikey"]
            pdb_df["match_key"] = pdb_df["pdb_inchikey"]

        # Get unique PDB ligand entries for matching
        pdb_unique = pdb_df.drop_duplicates(subset=["match_key", "pdb_id", "ligand_code"])

        # Merge activities with PDB ligands based on InChIKey
        linked = activities.merge(
            pdb_unique[["match_key", "pdb_id", "ligand_code", "pdb_inchikey"]],
            on="match_key",
            how="inner",
        )

        ligand_matched_count = len(linked)
        logger.info(f"InChIKey matching found {ligand_matched_count} potential links")

        # Validate that PDB structure contains BOTH protein AND ligand
        if has_protein_linking and len(linked) > 0:
            logger.info("Validating protein-ligand pairs (ensuring PDB contains both)")

            # Filter to keep only rows where pdb_id is in pdb_ids from protein linking
            linked["_valid_pair"] = linked.apply(
                lambda row: self._validate_protein_ligand_pair(row["pdb_id"], row["pdb_ids"]),
                axis=1,
            )

            validated_count = linked["_valid_pair"].sum()
            filtered_count = ligand_matched_count - validated_count

            logger.info(
                f"Validated {validated_count} protein-ligand pairs "
                f"(filtered out {filtered_count} where PDB doesn't contain both)"
            )

            # Keep only validated pairs
            linked = linked[linked["_valid_pair"]].drop(columns=["_valid_pair"])

        # Drop temporary columns
        linked = linked.drop(columns=["match_key"], errors="ignore")

        logger.info(
            f"Final result: {len(linked)} activities linked to PDB structures "
            f"containing both protein and ligand "
            f"(from {len(chembl_activities)} original activities)"
        )

        return linked

    def get_linkage_stats(
        self,
        chembl_compounds: pd.DataFrame,
        pdb_ligands: pd.DataFrame,
    ) -> dict:
        """Get statistics about ligand-level linkage.

        Args:
            chembl_compounds: ChEMBL compounds with InChIKeys (may also have pdb_ids
                             from protein-level linking)
            pdb_ligands: PDB ligands with InChIKeys

        Returns:
            Dictionary with linkage statistics including validated pair counts
        """
        # Find InChIKey columns
        chembl_inchikey_col = None
        for col in ["standard_inchi_key", "inchikey", "standard_inchikey"]:
            if col in chembl_compounds.columns:
                chembl_inchikey_col = col
                break

        pdb_inchikey_col = None
        for col in ["ligand_inchikey", "inchikey", "inchi_key"]:
            if col in pdb_ligands.columns:
                pdb_inchikey_col = col
                break

        if chembl_inchikey_col is None or pdb_inchikey_col is None:
            return {"error": "InChIKey columns not found"}

        chembl_keys = set(chembl_compounds[chembl_inchikey_col].dropna().unique())
        pdb_keys = set(pdb_ligands[pdb_inchikey_col].dropna().unique())

        if self.connectivity_only:
            chembl_keys = {k[:14] for k in chembl_keys if len(k) >= 14}
            pdb_keys = {k[:14] for k in pdb_keys if len(k) >= 14}

        common_keys = chembl_keys & pdb_keys

        stats = {
            "chembl_inchikey_count": len(chembl_keys),
            "pdb_inchikey_count": len(pdb_keys),
            "common_inchikey_count": len(common_keys),
            "match_type": ("connectivity_only" if self.connectivity_only else "full_inchikey"),
            "ligand_linkage_rate": (
                len(common_keys) / len(chembl_keys) if len(chembl_keys) > 0 else 0
            ),
        }

        # Add protein-ligand pair validation stats if protein linking info is available
        has_protein_linking = "pdb_ids" in chembl_compounds.columns
        stats["has_protein_linking"] = has_protein_linking

        if has_protein_linking and len(common_keys) > 0:
            # Count how many activities have at least one validated protein-ligand pair
            # by checking if any pdb_id from ligand matching is in pdb_ids from protein matching
            pdb_id_to_inchikeys: dict[str, set] = {}
            for _, row in pdb_ligands.iterrows():
                pdb_id = row.get("pdb_id")
                inchikey = row.get(pdb_inchikey_col)
                if pd.notna(pdb_id) and pd.notna(inchikey):
                    key = inchikey[:14] if self.connectivity_only else inchikey
                    if pdb_id not in pdb_id_to_inchikeys:
                        pdb_id_to_inchikeys[pdb_id] = set()
                    pdb_id_to_inchikeys[pdb_id].add(key)

            validated_activity_count = 0
            total_with_ligand_match = 0

            for _, row in chembl_compounds.iterrows():
                inchikey = row.get(chembl_inchikey_col)
                pdb_ids = row.get("pdb_ids")

                if pd.isna(inchikey) or pdb_ids is None:
                    continue

                key = inchikey[:14] if self.connectivity_only else inchikey
                if key not in common_keys:
                    continue

                total_with_ligand_match += 1

                # Check if any PDB from protein linking contains the ligand
                if hasattr(pdb_ids, "__iter__") and not isinstance(pdb_ids, str):
                    for pdb_id in pdb_ids:
                        pdb_id_upper = pdb_id.upper() if isinstance(pdb_id, str) else pdb_id
                        if pdb_id_upper in pdb_id_to_inchikeys:
                            if key in pdb_id_to_inchikeys[pdb_id_upper]:
                                validated_activity_count += 1
                                break

            stats["activities_with_ligand_match"] = total_with_ligand_match
            stats["activities_with_validated_pairs"] = validated_activity_count
            stats["validated_pair_rate"] = (
                validated_activity_count / total_with_ligand_match
                if total_with_ligand_match > 0
                else 0
            )

        return stats

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
