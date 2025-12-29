"""ChEMBL data downloader."""

import logging
import sqlite3
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional

import pandas as pd
from tqdm import tqdm

from chembl_pdb_linker.config import Config

logger = logging.getLogger(__name__)


class ChEMBLDownloader:
    """Download and extract ChEMBL data."""

    SQLITE_FILENAME_PATTERN = "chembl_*_sqlite.tar.gz"

    def __init__(self, config: Config):
        self.config = config
        self.paths = config.get_paths()
        self.raw_dir = self.paths["raw_dir"]
        self.intermediate_dir = self.paths["intermediate_dir"]

    def download_sqlite(self, version: str = "latest") -> Path:
        """Download ChEMBL SQLite database.

        Args:
            version: ChEMBL version (e.g., "34") or "latest"

        Returns:
            Path to the downloaded SQLite file
        """
        self.raw_dir.mkdir(parents=True, exist_ok=True)

        if version == "latest":
            ftp_url = f"{self.config.chembl.ftp_base}chembl_34_sqlite.tar.gz"
            filename = "chembl_34_sqlite.tar.gz"
        else:
            ftp_url = (
                f"ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/"
                f"releases/chembl_{version}/chembl_{version}_sqlite.tar.gz"
            )
            filename = f"chembl_{version}_sqlite.tar.gz"

        tar_path = self.raw_dir / filename
        sqlite_dir = self.raw_dir / f"chembl_{version}_sqlite"
        sqlite_path = sqlite_dir / f"chembl_{version}.db"

        if sqlite_path.exists():
            logger.info(f"ChEMBL SQLite already exists: {sqlite_path}")
            return sqlite_path

        if not tar_path.exists():
            logger.info(f"Downloading ChEMBL SQLite from {ftp_url}")
            self._download_with_progress(ftp_url, tar_path)

        logger.info(f"Extracting {tar_path}")
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=self.raw_dir)

        # Find the actual sqlite file
        for db_file in self.raw_dir.rglob("*.db"):
            if "chembl" in db_file.name.lower():
                return db_file

        raise FileNotFoundError("Could not find ChEMBL SQLite database after extraction")

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

    def extract_activities(self, sqlite_path: Path) -> pd.DataFrame:
        """Extract activity data from ChEMBL SQLite.

        Filters by:
        - Confidence score (default: 9)
        - Activity types (configurable)
        - Standard units

        Returns:
            DataFrame with activity data
        """
        logger.info("Extracting activities from ChEMBL")

        activity_types = self.config.chembl.activity_types
        confidence_score = self.config.chembl.confidence_score
        standard_units = self.config.chembl.standard_units

        activity_types_str = ",".join(f"'{t}'" for t in activity_types)
        units_str = ",".join(f"'{u}'" for u in standard_units)

        query = f"""
        SELECT
            act.activity_id,
            act.assay_id,
            act.molregno,
            act.standard_type,
            act.standard_value,
            act.standard_units,
            act.pchembl_value,
            mol.chembl_id AS molecule_chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi_key,
            ass.chembl_id AS assay_chembl_id,
            ass.confidence_score,
            td.chembl_id AS target_chembl_id,
            td.pref_name AS target_name,
            td.target_type,
            tc.accession AS uniprot_id
        FROM activities act
        JOIN molecule_dictionary mol ON act.molregno = mol.molregno
        JOIN compound_structures cs ON mol.molregno = cs.molregno
        JOIN assays ass ON act.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        LEFT JOIN target_components tc ON td.tid = tc.tid
        WHERE
            ass.confidence_score >= {confidence_score}
            AND act.standard_type IN ({activity_types_str})
            AND act.standard_units IN ({units_str})
            AND act.standard_value IS NOT NULL
            AND cs.standard_inchi_key IS NOT NULL
            AND tc.accession IS NOT NULL
        """

        conn = sqlite3.connect(sqlite_path)
        try:
            df = pd.read_sql_query(query, conn)
            logger.info(f"Extracted {len(df)} activity records")
            return df
        finally:
            conn.close()

    def extract_compounds(self, sqlite_path: Path) -> pd.DataFrame:
        """Extract compound structures with InChIKeys.

        Returns:
            DataFrame with compound data
        """
        logger.info("Extracting compounds from ChEMBL")

        query = """
        SELECT
            mol.molregno,
            mol.chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key
        FROM molecule_dictionary mol
        JOIN compound_structures cs ON mol.molregno = cs.molregno
        WHERE cs.standard_inchi_key IS NOT NULL
        """

        conn = sqlite3.connect(sqlite_path)
        try:
            df = pd.read_sql_query(query, conn)
            logger.info(f"Extracted {len(df)} compounds")
            return df
        finally:
            conn.close()

    def extract_targets(self, sqlite_path: Path) -> pd.DataFrame:
        """Extract target information with UniProt mappings.

        Returns:
            DataFrame with target data
        """
        logger.info("Extracting targets from ChEMBL")

        query = """
        SELECT DISTINCT
            td.tid,
            td.chembl_id AS target_chembl_id,
            td.pref_name AS target_name,
            td.target_type,
            td.organism,
            tc.accession AS uniprot_id
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        WHERE tc.accession IS NOT NULL
        """

        conn = sqlite3.connect(sqlite_path)
        try:
            df = pd.read_sql_query(query, conn)
            logger.info(f"Extracted {len(df)} target records")
            return df
        finally:
            conn.close()

    def save_intermediate(self, df: pd.DataFrame, name: str) -> Path:
        """Save intermediate data to parquet.

        Args:
            df: DataFrame to save
            name: Name for the output file

        Returns:
            Path to saved file
        """
        self.intermediate_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.intermediate_dir / f"chembl_{name}.parquet"
        df.to_parquet(output_path, compression="snappy")
        logger.info(f"Saved {name} to {output_path}")
        return output_path

    def run(self, version: str = "latest") -> dict[str, Path]:
        """Run the full ChEMBL download and extraction pipeline.

        Args:
            version: ChEMBL version to download

        Returns:
            Dictionary of output file paths
        """
        sqlite_path = self.download_sqlite(version)

        activities_df = self.extract_activities(sqlite_path)
        compounds_df = self.extract_compounds(sqlite_path)
        targets_df = self.extract_targets(sqlite_path)

        return {
            "activities": self.save_intermediate(activities_df, "activities"),
            "compounds": self.save_intermediate(compounds_df, "compounds"),
            "targets": self.save_intermediate(targets_df, "targets"),
        }
