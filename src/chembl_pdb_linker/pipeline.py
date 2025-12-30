"""Main pipeline orchestration for ChEMBL-PDB linking."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from chembl_pdb_linker.config import Config
from chembl_pdb_linker.downloaders.chembl import ChEMBLDownloader
from chembl_pdb_linker.downloaders.pdb import PDBDownloader
from chembl_pdb_linker.extractors.bioactivity import BioactivityExtractor
from chembl_pdb_linker.extractors.structure import StructureExtractor
from chembl_pdb_linker.linkers.ligand import LigandLinker
from chembl_pdb_linker.linkers.uniprot import UniProtLinker

logger = logging.getLogger(__name__)


class Pipeline:
    """Main pipeline for linking ChEMBL and PDB data."""

    def __init__(self, config: Optional[Config] = None, base_dir: Optional[Path] = None):
        """Initialize the pipeline.

        Args:
            config: Configuration object. If None, loads default config.
            base_dir: Base directory for the project.
        """
        if config is None:
            config = Config.default(base_dir)

        self.config = config
        self.paths = config.get_paths()

        # Initialize components
        self.chembl_downloader = ChEMBLDownloader(config)
        self.pdb_downloader = PDBDownloader(config)
        self.uniprot_linker = UniProtLinker(config)
        self.ligand_linker = LigandLinker(config)
        self.bioactivity_extractor = BioactivityExtractor(config)
        self.structure_extractor = StructureExtractor(config)

        # Ensure directories exist
        config.ensure_directories()

    def download(self, chembl_version: str = "latest") -> dict[str, Path]:
        """Download all required data.

        This coordinates ChEMBL and PDB data downloads to enable efficient
        ligand-structure linking via RCSB Search API.

        Args:
            chembl_version: ChEMBL version to download

        Returns:
            Dictionary of downloaded file paths
        """
        logger.info("=== Starting data download ===")

        outputs: dict[str, Path] = {}
        intermediate_dir = self.paths["intermediate_dir"]

        # Step 1: Download ChEMBL data first (needed for filtering PDB queries)
        logger.info("Step 1: Downloading ChEMBL data...")
        chembl_outputs = self.chembl_downloader.run(version=chembl_version)
        outputs.update({f"chembl_{k}": v for k, v in chembl_outputs.items()})

        # Step 2: Load ChEMBL activities to get InChIKeys and UniProt IDs
        logger.info("Step 2: Loading ChEMBL data for PDB filtering...")
        activities_path = intermediate_dir / "chembl_activities.parquet"
        if activities_path.exists():
            activities = pd.read_parquet(activities_path)
            chembl_inchikeys = set(activities["standard_inchi_key"].dropna().unique())
            chembl_uniprots = list(activities["uniprot_id"].dropna().unique())
            logger.info(
                f"ChEMBL data: {len(chembl_inchikeys):,} unique InChIKeys, "
                f"{len(chembl_uniprots):,} unique UniProt IDs"
            )
        else:
            logger.warning("ChEMBL activities not found, PDB download will be unfiltered")
            chembl_inchikeys = None
            chembl_uniprots = None

        # Step 3: Load or download ligand InChIKey mapping
        logger.info("Step 3: Loading ligand InChIKey mapping...")
        ligand_inchikey_mapping = self.pdb_downloader.download_ligand_inchikey_mapping()
        if ligand_inchikey_mapping.empty:
            logger.warning(
                "No ligand InChIKey mapping available. " "Ligand-level linking will be limited."
            )
            ligand_inchikey_mapping = None

        # Step 4: Download PDB/SIFTS data with efficient ligand lookup
        logger.info("Step 4: Downloading PDB/SIFTS data...")
        pdb_outputs = self.pdb_downloader.run(
            uniprot_ids=chembl_uniprots,
            chembl_inchikeys=chembl_inchikeys,
            ligand_inchikey_mapping=ligand_inchikey_mapping,
        )
        outputs.update({f"pdb_{k}": v for k, v in pdb_outputs.items()})

        logger.info(f"Download complete. Files: {list(outputs.keys())}")
        return outputs

    def link(
        self,
        activities_path: Optional[Path] = None,
        sifts_path: Optional[Path] = None,
        ligands_path: Optional[Path] = None,
    ) -> pd.DataFrame:
        """Link ChEMBL and PDB data.

        Args:
            activities_path: Path to ChEMBL activities parquet
            sifts_path: Path to SIFTS mapping parquet
            ligands_path: Path to PDB ligands parquet (optional)

        Returns:
            Linked DataFrame
        """
        logger.info("=== Starting data linking ===")

        intermediate_dir = self.paths["intermediate_dir"]

        # Load activities
        if activities_path is None:
            activities_path = intermediate_dir / "chembl_activities.parquet"
        if not activities_path.exists():
            raise FileNotFoundError(
                f"Activities file not found: {activities_path}. Run download first."
            )
        activities = pd.read_parquet(activities_path)
        logger.info(f"Loaded {len(activities)} activities")

        # Load SIFTS mapping
        if sifts_path is None:
            sifts_path = intermediate_dir / "pdb_sifts_mapping.parquet"
        if not sifts_path.exists():
            raise FileNotFoundError(f"SIFTS mapping not found: {sifts_path}. Run download first.")
        sifts = pd.read_parquet(sifts_path)
        logger.info(f"Loaded {len(sifts)} SIFTS mappings")

        # Step 1: Protein-level linking via UniProt
        logger.info("Step 1: Protein-level linking via UniProt")
        protein_stats = self.uniprot_linker.get_linkage_stats(activities, sifts)
        logger.info(f"Protein linkage stats: {protein_stats}")

        protein_linked = self.uniprot_linker.link(activities, sifts)
        logger.info(f"Protein-linked: {len(protein_linked)} activities")

        # Step 2: Ligand-level linking via InChIKey (if ligands available)
        if ligands_path is None:
            ligands_path = intermediate_dir / "pdb_ligands.parquet"

        if ligands_path.exists():
            logger.info("Step 2: Ligand-level linking via InChIKey")
            ligands = pd.read_parquet(ligands_path)
            logger.info(f"Loaded {len(ligands)} PDB ligands")

            ligand_stats = self.ligand_linker.get_linkage_stats(activities, ligands)
            logger.info(f"Ligand linkage stats: {ligand_stats}")

            # Link activities to PDB ligands
            ligand_linked = self.ligand_linker.link_with_activities(protein_linked, ligands)
            logger.info(f"Ligand-linked: {len(ligand_linked)} activities")

            # Use ligand-linked if we got matches, otherwise fall back to protein-linked
            if len(ligand_linked) > 0:
                linked = ligand_linked
            else:
                logger.warning("No ligand matches found, using protein-level linking only")
                linked = protein_linked
        else:
            logger.info("No PDB ligands file found, using protein-level linking only")
            linked = protein_linked

        # Save intermediate result
        self.uniprot_linker.save_intermediate(linked, "chembl_pdb")

        return linked

    def extract(
        self,
        linked_path: Optional[Path] = None,
        enrich_structures: bool = True,
    ) -> Path:
        """Extract and format the final dataset.

        Args:
            linked_path: Path to linked data parquet
            enrich_structures: Whether to fetch additional structure metadata

        Returns:
            Path to the final output file
        """
        logger.info("=== Starting data extraction ===")

        intermediate_dir = self.paths["intermediate_dir"]

        # Load linked data
        if linked_path is None:
            linked_path = intermediate_dir / "linked_chembl_pdb.parquet"
        if not linked_path.exists():
            raise FileNotFoundError(f"Linked data not found: {linked_path}. Run link first.")

        linked = pd.read_parquet(linked_path)
        logger.info(f"Loaded {len(linked)} linked records")

        # Standardize bioactivity values
        logger.info("Standardizing bioactivity values")
        linked = self.bioactivity_extractor.standardize_units(linked)

        # Compute pChEMBL if not present
        if "pchembl_value" not in linked.columns or bool(linked["pchembl_value"].isna().all()):
            logger.info("Computing pChEMBL values")
            linked = self.bioactivity_extractor.compute_pchembl(linked)

        # Enrich with structure metadata
        if enrich_structures:
            logger.info("Enriching with structure metadata")
            linked = self.structure_extractor.enrich_linked_data(linked, fetch_metadata=True)

        # Filter by resolution
        if "resolution" in linked.columns:
            linked = self.structure_extractor.filter_by_resolution(linked)

        # Extract final columns
        final_df = self.bioactivity_extractor.extract_for_linked_data(linked)

        # Save output
        output_path = self.bioactivity_extractor.save_output(final_df)

        logger.info("=== Pipeline complete ===")
        logger.info(f"Final dataset: {len(final_df)} records")
        logger.info(f"Output file: {output_path}")

        return output_path

    def run(
        self,
        chembl_version: str = "latest",
        enrich_structures: bool = True,
    ) -> Path:
        """Run the complete pipeline.

        Args:
            chembl_version: ChEMBL version to download
            enrich_structures: Whether to fetch additional structure metadata

        Returns:
            Path to the final output file
        """
        logger.info("=== Running complete ChEMBL-PDB linking pipeline ===")

        # Step 1: Download
        self.download(chembl_version=chembl_version)

        # Step 2: Link
        self.link()

        # Step 3: Extract
        output_path = self.extract(enrich_structures=enrich_structures)

        return output_path

    def get_statistics(self) -> dict:
        """Get statistics about the current dataset.

        Returns:
            Dictionary with dataset statistics
        """
        intermediate_dir = self.paths["intermediate_dir"]
        output_dir = self.paths["output_dir"]

        stats = {}

        # Check intermediate files
        activities_path = intermediate_dir / "chembl_activities.parquet"
        if activities_path.exists():
            df = pd.read_parquet(activities_path)
            stats["chembl_activities"] = len(df)
            stats["chembl_unique_compounds"] = df["molecule_chembl_id"].nunique()
            stats["chembl_unique_targets"] = df["target_chembl_id"].nunique()

        sifts_path = intermediate_dir / "pdb_sifts_mapping.parquet"
        if sifts_path.exists():
            df = pd.read_parquet(sifts_path)
            stats["pdb_structures"] = df["pdb_id"].nunique()
            stats["pdb_uniprot_ids"] = df["uniprot_id"].nunique()

        linked_path = intermediate_dir / "linked_chembl_pdb.parquet"
        if linked_path.exists():
            df = pd.read_parquet(linked_path)
            stats["linked_records"] = len(df)

        # Check output files
        output_path = output_dir / "bioactivity_pdb_linked.parquet"
        if output_path.exists():
            df = pd.read_parquet(output_path)
            stats["final_records"] = len(df)
            if "chembl_id" in df.columns:
                stats["final_unique_compounds"] = df["chembl_id"].nunique()
            if "pdb_id" in df.columns:
                stats["final_unique_structures"] = df["pdb_id"].nunique()

        return stats
