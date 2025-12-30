"""Configuration handling for ChEMBL-PDB Linker."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass
class PathsConfig:
    """Paths configuration."""

    data_dir: str = "data"
    raw_dir: str = "data/raw"
    intermediate_dir: str = "data/intermediate"
    output_dir: str = "data/curated"

    def as_paths(self, base_dir: Path) -> dict[str, Path]:
        """Convert to absolute Path objects."""
        return {
            "data_dir": base_dir / self.data_dir,
            "raw_dir": base_dir / self.raw_dir,
            "intermediate_dir": base_dir / self.intermediate_dir,
            "output_dir": base_dir / self.output_dir,
        }


@dataclass
class ChEMBLConfig:
    """ChEMBL-specific configuration."""

    ftp_base: str = "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    api_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    confidence_score: int = 9
    activity_types: list[str] = field(default_factory=lambda: ["IC50", "Ki", "Kd", "EC50"])
    standard_units: list[str] = field(default_factory=lambda: ["nM"])


@dataclass
class PDBConfig:
    """PDB-specific configuration."""

    api_base: str = "https://www.ebi.ac.uk/pdbe/api"
    rcsb_api_base: str = "https://data.rcsb.org/rest/v1"
    rcsb_search_api: str = "https://search.rcsb.org/rcsbsearch/v2/query"
    sifts_ftp: str = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"
    max_resolution: float = 3.5


@dataclass
class LinkingConfig:
    """Linking-specific configuration."""

    use_inchikey: bool = True
    inchikey_connectivity_only: bool = False


@dataclass
class APIConfig:
    """API settings configuration."""

    rate_limit: int = 5
    timeout: int = 30
    max_retries: int = 3


@dataclass
class OutputConfig:
    """Output settings configuration."""

    format: str = "parquet"
    compression: str = "snappy"


@dataclass
class Config:
    """Main configuration class."""

    paths: PathsConfig = field(default_factory=PathsConfig)
    chembl: ChEMBLConfig = field(default_factory=ChEMBLConfig)
    pdb: PDBConfig = field(default_factory=PDBConfig)
    linking: LinkingConfig = field(default_factory=LinkingConfig)
    api: APIConfig = field(default_factory=APIConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    base_dir: Path = field(default_factory=Path.cwd)

    @classmethod
    def from_yaml(cls, config_path: Path, base_dir: Optional[Path] = None) -> "Config":
        """Load configuration from a YAML file."""
        with open(config_path) as f:
            data = yaml.safe_load(f)

        if base_dir is None:
            base_dir = config_path.parent.parent

        config = cls(base_dir=base_dir)

        if "paths" in data:
            config.paths = PathsConfig(**data["paths"])
        if "chembl" in data:
            config.chembl = ChEMBLConfig(**data["chembl"])
        if "pdb" in data:
            config.pdb = PDBConfig(**data["pdb"])
        if "linking" in data:
            config.linking = LinkingConfig(**data["linking"])
        if "api" in data:
            config.api = APIConfig(**data["api"])
        if "output" in data:
            config.output = OutputConfig(**data["output"])

        return config

    @classmethod
    def default(cls, base_dir: Optional[Path] = None) -> "Config":
        """Load the default configuration."""
        if base_dir is None:
            base_dir = Path.cwd()

        default_config = base_dir / "config" / "default.yaml"
        if default_config.exists():
            return cls.from_yaml(default_config, base_dir)
        return cls(base_dir=base_dir)

    def get_paths(self) -> dict[str, Path]:
        """Get all paths as absolute Path objects."""
        return self.paths.as_paths(self.base_dir)

    def ensure_directories(self) -> None:
        """Create all necessary directories."""
        for path in self.get_paths().values():
            path.mkdir(parents=True, exist_ok=True)
