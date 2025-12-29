"""PDB data downloader."""

import gzip
import logging
import time
import urllib.request
from io import StringIO
from pathlib import Path
from typing import Optional

import httpx
import pandas as pd
from tqdm import tqdm

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class PDBDownloader:
    """Download and process PDB/PDBe data."""

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.raw_dir = self.paths["raw_dir"]
        self.intermediate_dir = self.paths["intermediate_dir"]
        self.client = httpx.Client(timeout=config.api.timeout)

    def __del__(self):
        if hasattr(self, "client"):
            self.client.close()

    def download_sifts_mapping(self) -> Path:
        """Download SIFTS UniProt-to-PDB mappings.

        Returns:
            Path to the downloaded mapping file
        """
        self.raw_dir.mkdir(parents=True, exist_ok=True)

        url = f"{self.config.pdb.sifts_ftp}pdb_chain_uniprot.csv.gz"
        gz_path = self.raw_dir / "pdb_chain_uniprot.csv.gz"
        csv_path = self.raw_dir / "pdb_chain_uniprot.csv"

        if csv_path.exists():
            logger.info(f"SIFTS mapping already exists: {csv_path}")
            return csv_path

        if not gz_path.exists():
            logger.info(f"Downloading SIFTS mapping from {url}")
            self._download_with_progress(url, gz_path)

        logger.info(f"Extracting {gz_path}")
        with gzip.open(gz_path, "rt") as f_in:
            with open(csv_path, "w") as f_out:
                f_out.write(f_in.read())

        return csv_path

    def _download_with_progress(self, url: str, dest: Path) -> None:
        """Download a file with progress bar."""

        class DownloadProgressBar(tqdm):
            def update_to(self, b: int = 1, bsize: int = 1, tsize: Optional[int] = None) -> None:
                if tsize is not None:
                    self.total = tsize
                self.update(b * bsize - self.n)

        with DownloadProgressBar(
            unit="B", unit_scale=True, miniters=1, desc=dest.name
        ) as progress:
            urllib.request.urlretrieve(url, dest, reporthook=progress.update_to)

    def load_sifts_mapping(self, csv_path: Optional[Path] = None) -> pd.DataFrame:
        """Load and parse SIFTS UniProt-to-PDB mappings.

        Returns:
            DataFrame with PDB-UniProt mappings
        """
        if csv_path is None:
            csv_path = self.download_sifts_mapping()

        logger.info(f"Loading SIFTS mapping from {csv_path}")

        # SIFTS CSV has a header line starting with #
        df = pd.read_csv(csv_path, comment="#", header=0)

        # Rename columns to be more descriptive
        column_mapping = {
            "PDB": "pdb_id",
            "CHAIN": "chain_id",
            "SP_PRIMARY": "uniprot_id",
            "RES_BEG": "res_begin",
            "RES_END": "res_end",
            "PDB_BEG": "pdb_begin",
            "PDB_END": "pdb_end",
            "SP_BEG": "sp_begin",
            "SP_END": "sp_end",
        }
        df = df.rename(columns=column_mapping)

        logger.info(f"Loaded {len(df)} PDB-UniProt mappings")
        return df

    def fetch_pdb_entry_info(self, pdb_ids: list[str]) -> pd.DataFrame:
        """Fetch metadata for PDB entries using PDBe API.

        Args:
            pdb_ids: List of PDB IDs to fetch

        Returns:
            DataFrame with PDB entry metadata
        """
        logger.info(f"Fetching metadata for {len(pdb_ids)} PDB entries")

        results = []
        rate_limit_delay = 1.0 / self.config.api.rate_limit

        for i, pdb_id in enumerate(tqdm(pdb_ids, desc="Fetching PDB metadata")):
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
                                "resolution": entry.get("resolution", None),
                                "experimental_method": entry.get(
                                    "experimental_method", [""]
                                )[0],
                            }
                        )
                elif response.status_code != 404:
                    logger.warning(f"Failed to fetch {pdb_id}: {response.status_code}")

            except Exception as e:
                logger.warning(f"Error fetching {pdb_id}: {e}")

            # Rate limiting
            if i < len(pdb_ids) - 1:
                time.sleep(rate_limit_delay)

        return pd.DataFrame(results)

    def fetch_ligand_info(self, pdb_ids: list[str]) -> pd.DataFrame:
        """Fetch ligand information for PDB entries.

        Args:
            pdb_ids: List of PDB IDs

        Returns:
            DataFrame with ligand information
        """
        logger.info(f"Fetching ligand info for {len(pdb_ids)} PDB entries")

        results = []
        rate_limit_delay = 1.0 / self.config.api.rate_limit

        for i, pdb_id in enumerate(tqdm(pdb_ids, desc="Fetching ligand info")):
            try:
                url = f"{self.config.pdb.api_base}/pdb/entry/ligand_monomers/{pdb_id.lower()}"
                response = self.client.get(url)

                if response.status_code == 200:
                    data = response.json()
                    if pdb_id.lower() in data:
                        for ligand in data[pdb_id.lower()]:
                            results.append(
                                {
                                    "pdb_id": pdb_id.upper(),
                                    "ligand_code": ligand.get("chem_comp_id", ""),
                                    "ligand_name": ligand.get("chem_comp_name", ""),
                                    "chain_id": ligand.get("chain_id", ""),
                                }
                            )

            except Exception as e:
                logger.warning(f"Error fetching ligands for {pdb_id}: {e}")

            if i < len(pdb_ids) - 1:
                time.sleep(rate_limit_delay)

        return pd.DataFrame(results)

    def download_ligand_expo(self) -> Path:
        """Download the PDB Chemical Component Dictionary (Ligand Expo).

        This contains InChIKeys for PDB ligands.

        Returns:
            Path to the downloaded file
        """
        self.raw_dir.mkdir(parents=True, exist_ok=True)

        url = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
        gz_path = self.raw_dir / "components.cif.gz"

        if gz_path.exists():
            logger.info(f"Component dictionary already exists: {gz_path}")
            return gz_path

        logger.info(f"Downloading PDB Component Dictionary from {url}")
        self._download_with_progress(url, gz_path)
        return gz_path

    def fetch_ligand_inchikeys_batch(self, ligand_codes: list[str]) -> pd.DataFrame:
        """Fetch InChIKeys for PDB ligands via PDBe API.

        Args:
            ligand_codes: List of 3-letter ligand codes

        Returns:
            DataFrame with ligand codes and InChIKeys
        """
        logger.info(f"Fetching InChIKeys for {len(ligand_codes)} ligands")

        results = []
        rate_limit_delay = 1.0 / self.config.api.rate_limit
        unique_codes = list(set(ligand_codes))

        for i, code in enumerate(tqdm(unique_codes, desc="Fetching ligand InChIKeys")):
            try:
                url = f"{self.config.pdb.api_base}/pdb/compound/summary/{code}"
                response = self.client.get(url)

                if response.status_code == 200:
                    data = response.json()
                    if code in data:
                        entry = data[code][0]
                        inchikey = entry.get("inchi_key", "")
                        if inchikey:
                            results.append(
                                {
                                    "ligand_code": code,
                                    "ligand_inchikey": inchikey,
                                    "ligand_name": entry.get("name", ""),
                                    "ligand_formula": entry.get("formula", ""),
                                }
                            )

            except Exception as e:
                logger.warning(f"Error fetching InChIKey for {code}: {e}")

            if i < len(unique_codes) - 1:
                time.sleep(rate_limit_delay)

        return pd.DataFrame(results)

    def save_intermediate(self, df: pd.DataFrame, name: str) -> Path:
        """Save intermediate data to parquet.

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

    def run(self, uniprot_ids: Optional[list[str]] = None) -> dict[str, Path]:
        """Run the PDB data download and extraction pipeline.

        Args:
            uniprot_ids: Optional list of UniProt IDs to filter by.
                        If None, downloads full SIFTS mapping.

        Returns:
            Dictionary of output file paths
        """
        # Download and load SIFTS mapping
        sifts_path = self.download_sifts_mapping()
        sifts_df = self.load_sifts_mapping(sifts_path)

        # Filter by UniProt IDs if provided
        if uniprot_ids:
            sifts_df = sifts_df[sifts_df["uniprot_id"].isin(uniprot_ids)]
            logger.info(f"Filtered to {len(sifts_df)} mappings for {len(uniprot_ids)} UniProt IDs")

        # Get unique PDB IDs
        pdb_ids = sifts_df["pdb_id"].unique().tolist()
        logger.info(f"Found {len(pdb_ids)} unique PDB IDs")

        # Fetch PDB metadata (for a subset if too many)
        if len(pdb_ids) > 10000:
            logger.warning(
                f"Too many PDB IDs ({len(pdb_ids)}). "
                "Consider filtering by UniProt IDs first."
            )

        # Save mappings
        outputs = {
            "sifts_mapping": self.save_intermediate(sifts_df, "sifts_mapping"),
        }

        # Optionally fetch additional metadata via API
        # This is expensive for large datasets, so we make it optional
        if len(pdb_ids) <= 1000:
            ligands_df = self.fetch_ligand_info(pdb_ids)
            if not ligands_df.empty:
                # Fetch InChIKeys for ligands
                ligand_codes = ligands_df["ligand_code"].unique().tolist()
                inchikeys_df = self.fetch_ligand_inchikeys_batch(ligand_codes)

                # Merge InChIKeys with ligand info
                ligands_df = ligands_df.merge(inchikeys_df, on="ligand_code", how="left")
                outputs["ligands"] = self.save_intermediate(ligands_df, "ligands")

        return outputs
