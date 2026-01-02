"""Structure metadata extraction from PDB."""

import logging
import time
from pathlib import Path
from typing import Optional

import httpx
import pandas as pd
from tqdm import tqdm

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class StructureExtractor:
    """Extract and process structure metadata from PDB."""

    # Base URLs for structure downloads
    RCSB_DOWNLOAD_BASE = "https://files.rcsb.org/download"
    PDBE_DOWNLOAD_BASE = "https://www.ebi.ac.uk/pdbe/entry-files/download"

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.intermediate_dir = self.paths["intermediate_dir"]
        self.output_dir = self.paths["output_dir"]
        self.client = httpx.Client(timeout=config.api.timeout)

    def __del__(self):
        if hasattr(self, "client"):
            self.client.close()

    def fetch_entry_metadata(self, pdb_ids: list[str]) -> pd.DataFrame:
        """Fetch metadata for PDB entries.

        Args:
            pdb_ids: List of PDB IDs

        Returns:
            DataFrame with structure metadata
        """
        logger.info(f"Fetching metadata for {len(pdb_ids)} PDB entries")

        results = []
        rate_limit_delay = 1.0 / self.config.api.rate_limit
        unique_ids = list(set(pdb_ids))

        for i, pdb_id in enumerate(tqdm(unique_ids, desc="Fetching PDB metadata")):
            try:
                url = f"{self.config.pdb.api_base}/pdb/entry/summary/{pdb_id.lower()}"
                response = self.client.get(url)

                if response.status_code == 200:
                    data = response.json()
                    if pdb_id.lower() in data:
                        entry = data[pdb_id.lower()][0]
                        results.append(
                            {
                                "pdb_id": pdb_id.upper(),
                                "title": entry.get("title", ""),
                                "release_date": entry.get("release_date", ""),
                                "deposition_date": entry.get("deposition_date", ""),
                                "revision_date": entry.get("revision_date", ""),
                                "resolution": entry.get("resolution", None),
                                "experimental_method": (
                                    entry.get("experimental_method", [""])[0]
                                    if entry.get("experimental_method")
                                    else ""
                                ),
                                "number_of_entities": entry.get("number_of_entities", {}).get(
                                    "polypeptide(L)", 0
                                ),
                            }
                        )

            except Exception as e:
                logger.warning(f"Error fetching metadata for {pdb_id}: {e}")

            if i < len(unique_ids) - 1:
                time.sleep(rate_limit_delay)

        return pd.DataFrame(results)

    def fetch_binding_sites(self, pdb_ids: list[str]) -> pd.DataFrame:
        """Fetch binding site information for PDB entries.

        Args:
            pdb_ids: List of PDB IDs

        Returns:
            DataFrame with binding site information
        """
        logger.info(f"Fetching binding sites for {len(pdb_ids)} PDB entries")

        results = []
        rate_limit_delay = 1.0 / self.config.api.rate_limit
        unique_ids = list(set(pdb_ids))

        for i, pdb_id in enumerate(tqdm(unique_ids, desc="Fetching binding sites")):
            try:
                url = f"{self.config.pdb.api_base}/pdb/entry/binding_sites/{pdb_id.lower()}"
                response = self.client.get(url)

                if response.status_code == 200:
                    data = response.json()
                    if pdb_id.lower() in data:
                        for site in data[pdb_id.lower()]:
                            results.append(
                                {
                                    "pdb_id": pdb_id.upper(),
                                    "site_id": site.get("site_id", ""),
                                    "ligand_code": site.get("ligand_residue_name", ""),
                                    "chain_id": site.get("chain_id", ""),
                                    "residue_count": len(site.get("site_residues", [])),
                                }
                            )

            except Exception as e:
                logger.debug(f"No binding sites for {pdb_id}: {e}")

            if i < len(unique_ids) - 1:
                time.sleep(rate_limit_delay)

        return pd.DataFrame(results)

    def generate_download_urls(
        self,
        pdb_ids: list[str],
        file_format: str = "cif",
    ) -> pd.DataFrame:
        """Generate download URLs for structure files.

        Args:
            pdb_ids: List of PDB IDs
            file_format: File format (cif, pdb, pdb.gz, cif.gz)

        Returns:
            DataFrame with PDB IDs and download URLs
        """
        logger.info(f"Generating download URLs for {len(pdb_ids)} structures")

        results = []
        for pdb_id in pdb_ids:
            pdb_lower = pdb_id.lower()

            # RCSB URL
            if file_format == "cif":
                rcsb_url = f"{self.RCSB_DOWNLOAD_BASE}/{pdb_lower}.cif"
            elif file_format == "pdb":
                rcsb_url = f"{self.RCSB_DOWNLOAD_BASE}/{pdb_lower}.pdb"
            else:
                rcsb_url = f"{self.RCSB_DOWNLOAD_BASE}/{pdb_lower}.{file_format}"

            # PDBe URL
            pdbe_url = f"{self.PDBE_DOWNLOAD_BASE}/{pdb_lower}.cif"

            results.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "rcsb_url": rcsb_url,
                    "pdbe_url": pdbe_url,
                }
            )

        return pd.DataFrame(results)

    def filter_by_resolution(
        self,
        metadata: pd.DataFrame,
        max_resolution: Optional[float] = None,
    ) -> pd.DataFrame:
        """Filter structures by resolution.

        Args:
            metadata: DataFrame with resolution column
            max_resolution: Maximum resolution in Angstroms

        Returns:
            Filtered DataFrame
        """
        if max_resolution is None:
            max_resolution = self.config.pdb.max_resolution

        if "resolution" not in metadata.columns:
            logger.warning("No resolution column found")
            return metadata

        df = metadata.copy()
        initial_count = len(df)

        # Keep entries with resolution <= max_resolution or null resolution (NMR, etc.)
        df = df[(df["resolution"].isna()) | (df["resolution"] <= max_resolution)]

        logger.info(
            f"Filtered from {initial_count} to {len(df)} entries (resolution <= {max_resolution} Å)"
        )

        return df  # type: ignore[return-value]

    def enrich_linked_data(
        self,
        linked_data: pd.DataFrame,
        fetch_metadata: bool = True,
    ) -> pd.DataFrame:
        """Enrich linked data with structure metadata.

        Args:
            linked_data: DataFrame with linked ChEMBL-PDB data
            fetch_metadata: Whether to fetch additional metadata via API

        Returns:
            Enriched DataFrame
        """
        logger.info("Enriching linked data with structure metadata")

        df = linked_data.copy()

        # Extract PDB IDs
        pdb_ids = []
        if "pdb_id" in df.columns:
            pdb_ids = df["pdb_id"].dropna().unique().tolist()
        elif "pdb_ids" in df.columns:
            # Handle list/array column
            for ids in df["pdb_ids"].dropna():
                if hasattr(ids, "__iter__") and not isinstance(ids, str):
                    pdb_ids.extend(list(ids))
            pdb_ids = list(set(pdb_ids))

        if not pdb_ids:
            logger.warning("No PDB IDs found in linked data")
            return df

        logger.info(f"Found {len(pdb_ids)} unique PDB IDs")

        # Generate download URLs
        urls_df = self.generate_download_urls(pdb_ids)

        # Optionally fetch metadata
        if fetch_metadata and len(pdb_ids) <= 1000:
            metadata_df = self.fetch_entry_metadata(pdb_ids)

            # Merge URLs with metadata
            urls_df = urls_df.merge(metadata_df, on="pdb_id", how="left")

        # Merge with linked data
        if "pdb_id" in df.columns:
            df = df.merge(urls_df, on="pdb_id", how="left")
        else:
            # For data with pdb_ids as lists/arrays, we'll just add URLs for the first PDB
            df["primary_pdb_id"] = df["pdb_ids"].apply(
                lambda x: (
                    x[0]
                    if hasattr(x, "__len__") and not isinstance(x, str) and len(x) > 0
                    else None
                )
            )
            urls_df = urls_df.rename(
                columns={
                    "pdb_id": "primary_pdb_id",
                    "rcsb_url": "primary_rcsb_url",
                    "pdbe_url": "primary_pdbe_url",
                }
            )
            df = df.merge(urls_df, on="primary_pdb_id", how="left")

        logger.info(f"Enriched {len(df)} records with structure metadata")

        return df

    def save_intermediate(self, df: pd.DataFrame, name: str) -> Path:
        """Save intermediate data.

        Args:
            df: DataFrame to save
            name: Name for the output file

        Returns:
            Path to saved file
        """
        self.intermediate_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.intermediate_dir / f"pdb_{name}.parquet"
        df.to_parquet(output_path, compression="snappy")
        logger.info(f"Saved {name} to {output_path}")
        return output_path
