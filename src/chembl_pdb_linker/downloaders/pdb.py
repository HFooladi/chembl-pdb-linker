"""PDB data downloader."""

import gzip
import logging
import time
import urllib.request
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

        with DownloadProgressBar(unit="B", unit_scale=True, miniters=1, desc=dest.name) as progress:
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
        df = pd.read_csv(csv_path, comment="#", header=0, low_memory=False)

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

        # Convert string columns to string type and numeric columns properly
        str_cols = ["pdb_id", "chain_id", "uniprot_id"]
        for col in str_cols:
            if col in df.columns:
                df[col] = df[col].astype(str)

        # Convert numeric columns, coercing errors to NaN
        num_cols = ["res_begin", "res_end", "pdb_begin", "pdb_end", "sp_begin", "sp_end"]
        for col in num_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

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
                                "experimental_method": entry.get("experimental_method", [""])[0],
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

    def fetch_ligand_structures_via_rcsb(
        self,
        ligand_codes: list[str],
        filter_pdb_ids: Optional[set[str]] = None,
    ) -> pd.DataFrame:
        """Fetch PDB structures containing specific ligands via RCSB Search API.

        This method efficiently queries which PDB structures contain specific ligands
        by using the RCSB Search API (reverse lookup by ligand code).

        Args:
            ligand_codes: List of PDB ligand codes to search for
            filter_pdb_ids: Optional set of PDB IDs to filter results to
                           (e.g., PDBs containing ChEMBL proteins)

        Returns:
            DataFrame with (pdb_id, ligand_code) pairs
        """
        logger.info(
            f"Fetching PDB structures for {len(ligand_codes)} ligand codes via RCSB Search API"
        )

        results = []
        search_url = self.config.pdb.rcsb_search_api

        for code in tqdm(ligand_codes, desc="Querying RCSB for ligand structures"):
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text_chem",
                    "parameters": {
                        "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                        "operator": "exact_match",
                        "value": code,
                    },
                },
                "return_type": "entry",
                "request_options": {"return_all_hits": True},
            }

            try:
                response = self.client.post(search_url, json=query)
                if response.status_code == 200:
                    data = response.json()
                    for result in data.get("result_set", []):
                        pdb_id = result["identifier"].upper()
                        # Filter to relevant PDBs if filter set provided
                        if filter_pdb_ids is None or pdb_id in filter_pdb_ids:
                            results.append({"pdb_id": pdb_id, "ligand_code": code})
                time.sleep(0.05)  # Rate limiting
            except Exception as e:
                logger.debug(f"Error querying RCSB for {code}: {e}")

        df = pd.DataFrame(results)
        if not df.empty:
            df = df.drop_duplicates()
            logger.info(
                f"Found {len(df)} PDB-ligand pairs "
                f"({df['pdb_id'].nunique()} unique PDBs, "
                f"{df['ligand_code'].nunique()} unique ligands)"
            )
        else:
            logger.warning("No PDB-ligand pairs found")

        return df

    def fetch_pdb_ligands_with_inchikeys(
        self,
        chembl_inchikeys: set[str],
        chembl_pdb_ids: set[str],
        ligand_inchikey_mapping: pd.DataFrame,
    ) -> pd.DataFrame:
        """Fetch complete PDB ligand data with InChIKeys for linking.

        This method:
        1. Finds ligand codes that have InChIKeys matching ChEMBL compounds
        2. Queries RCSB to find which PDB structures contain those ligands
        3. Filters to PDB structures that contain ChEMBL proteins
        4. Returns complete mapping with InChIKeys for linking

        Args:
            chembl_inchikeys: Set of InChIKeys from ChEMBL compounds
            chembl_pdb_ids: Set of PDB IDs that contain ChEMBL proteins
            ligand_inchikey_mapping: DataFrame with ligand_code -> inchikey mapping

        Returns:
            DataFrame with (pdb_id, ligand_code, ligand_inchikey) for linking
        """
        logger.info("Building PDB ligand data for ChEMBL-PDB linking")

        # Find common InChIKeys between ChEMBL and PDB ligands
        inchikey_col = None
        for col in ["ligand_inchikey", "pdb_inchikey", "inchikey"]:
            if col in ligand_inchikey_mapping.columns:
                inchikey_col = col
                break

        if inchikey_col is None:
            raise ValueError("No InChIKey column found in ligand mapping")

        ligand_code_col = None
        for col in ["ligand_code", "pdb_ligand_code"]:
            if col in ligand_inchikey_mapping.columns:
                ligand_code_col = col
                break

        if ligand_code_col is None:
            raise ValueError("No ligand code column found in ligand mapping")

        # Find ligand codes with InChIKeys that match ChEMBL compounds
        pdb_inchikeys = set(ligand_inchikey_mapping[inchikey_col].dropna().unique())
        common_inchikeys = chembl_inchikeys & pdb_inchikeys
        logger.info(f"Common InChIKeys between ChEMBL and PDB: {len(common_inchikeys)}")

        if len(common_inchikeys) == 0:
            logger.warning("No common InChIKeys found")
            return pd.DataFrame()

        # Get ligand codes for common InChIKeys
        common_ligands = ligand_inchikey_mapping[
            ligand_inchikey_mapping[inchikey_col].isin(common_inchikeys)
        ]
        ligand_codes = common_ligands[ligand_code_col].unique().tolist()
        logger.info(f"Ligand codes to query: {len(ligand_codes)}")

        # Normalize PDB IDs to uppercase for matching
        chembl_pdb_ids_upper = {p.upper() for p in chembl_pdb_ids}

        # Query RCSB for PDB structures containing these ligands
        pdb_ligand_pairs = self.fetch_ligand_structures_via_rcsb(
            ligand_codes, filter_pdb_ids=chembl_pdb_ids_upper
        )

        if pdb_ligand_pairs.empty:
            return pd.DataFrame()

        # Add InChIKeys to the results
        ligand_to_inchikey = common_ligands[[ligand_code_col, inchikey_col]].drop_duplicates()
        ligand_to_inchikey = ligand_to_inchikey.rename(
            columns={ligand_code_col: "ligand_code", inchikey_col: "ligand_inchikey"}
        )

        result = pdb_ligand_pairs.merge(ligand_to_inchikey, on="ligand_code", how="left")

        logger.info(
            f"Final PDB ligand data: {len(result)} pairs "
            f"({result['pdb_id'].nunique()} PDBs, {result['ligand_code'].nunique()} ligands)"
        )

        return result

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

    def run(
        self,
        uniprot_ids: Optional[list[str]] = None,
        chembl_inchikeys: Optional[set[str]] = None,
        ligand_inchikey_mapping: Optional[pd.DataFrame] = None,
    ) -> dict[str, Path]:
        """Run the PDB data download and extraction pipeline.

        Args:
            uniprot_ids: Optional list of UniProt IDs to filter by.
                        If None, downloads full SIFTS mapping.
            chembl_inchikeys: Optional set of InChIKeys from ChEMBL compounds.
                             Used for efficient ligand-structure lookup via RCSB.
            ligand_inchikey_mapping: Optional DataFrame with ligand code to InChIKey mapping.
                                    If not provided, will attempt to load from raw data.

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
        pdb_ids = list(sifts_df["pdb_id"].unique())  # type: ignore[attr-defined]
        pdb_ids_set = {p.upper() for p in pdb_ids}
        logger.info(f"Found {len(pdb_ids)} unique PDB IDs")

        # Save mappings
        outputs: dict[str, Path] = {
            "sifts_mapping": self.save_intermediate(sifts_df, "sifts_mapping"),  # type: ignore[arg-type]
        }

        # Fetch PDB ligand data for linking
        if chembl_inchikeys is not None and ligand_inchikey_mapping is not None:
            # Use RCSB Search API for efficient ligand-structure lookup
            logger.info("Fetching PDB ligand data via RCSB Search API")
            ligands_df = self.fetch_pdb_ligands_with_inchikeys(
                chembl_inchikeys=chembl_inchikeys,
                chembl_pdb_ids=pdb_ids_set,
                ligand_inchikey_mapping=ligand_inchikey_mapping,
            )
            if not ligands_df.empty:
                outputs["ligands"] = self.save_intermediate(ligands_df, "ligands")
        elif len(pdb_ids) <= 1000:
            # Fallback: fetch ligand info per PDB ID for small datasets
            logger.info("Fetching PDB ligand info via PDBe API (small dataset)")
            ligands_df = self.fetch_ligand_info(pdb_ids)
            if not ligands_df.empty:
                # Fetch InChIKeys for ligands
                ligand_codes = ligands_df["ligand_code"].unique().tolist()
                inchikeys_df = self.fetch_ligand_inchikeys_batch(ligand_codes)

                # Merge InChIKeys with ligand info
                ligands_df = ligands_df.merge(inchikeys_df, on="ligand_code", how="left")
                outputs["ligands"] = self.save_intermediate(ligands_df, "ligands")
        else:
            logger.warning(
                f"Large dataset ({len(pdb_ids)} PDB IDs) but no ChEMBL InChIKeys provided. "
                "Ligand linking will not be available. "
                "Provide chembl_inchikeys and ligand_inchikey_mapping for full linking."
            )

        return outputs

    def download_ligand_inchikey_mapping(self) -> pd.DataFrame:
        """Download and parse ligand code to InChIKey mapping.

        This fetches InChIKeys for all common PDB ligands. For efficiency,
        it uses the PDBe compound summary API.

        Returns:
            DataFrame with (ligand_code, ligand_inchikey) mapping
        """
        # Check if we already have the mapping cached
        cache_path = self.intermediate_dir / "pdb_ligand_inchikeys.parquet"
        if cache_path.exists():
            logger.info(f"Loading cached ligand InChIKey mapping from {cache_path}")
            return pd.read_parquet(cache_path)

        # Also check raw directory for pre-downloaded mapping
        raw_cache = self.raw_dir / "pdb_ligand_inchikeys.parquet"
        if raw_cache.exists():
            logger.info(f"Loading ligand InChIKey mapping from {raw_cache}")
            df = pd.read_parquet(raw_cache)
            # Normalize column names
            if "pdb_ligand_code" in df.columns and "ligand_code" not in df.columns:
                df = df.rename(columns={"pdb_ligand_code": "ligand_code"})
            if "pdb_inchikey" in df.columns and "ligand_inchikey" not in df.columns:
                df = df.rename(columns={"pdb_inchikey": "ligand_inchikey"})
            return df

        logger.warning(
            "No ligand InChIKey mapping found. "
            "This mapping should be pre-generated for large-scale processing. "
            "Place pdb_ligand_inchikeys.parquet in data/raw/ directory."
        )
        return pd.DataFrame()
