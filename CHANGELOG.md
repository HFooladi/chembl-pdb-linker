# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-01-02

### Added

- Initial release of ChEMBL-PDB Linker
- **Download phase**: Automated download of ChEMBL SQLite database, SIFTS UniProt-PDB mappings, and PDB ligand data
- **Linking phase**:
  - Protein-level linking via UniProt IDs
  - Ligand-level linking via InChIKey with protein-ligand pair validation
  - RCSB Search API integration for efficient ligand-to-structure mapping
- **Extract phase**: Bioactivity standardization and structure metadata enrichment
- CLI interface with commands: `run`, `download`, `link`, `extract`, `stats`
- Configurable pipeline via YAML configuration files
- Support for multiple activity types: IC50, Ki, Kd, EC50
- Parquet output format for efficient data storage
- Comprehensive test suite with unit and integration tests
- CI/CD pipeline with GitHub Actions

### Output

- ~98,500 validated protein-ligand pairs with bioactivity data
- ~9,000 unique compounds
- ~14,700 unique PDB structures
- ~1,300 unique target proteins

[0.1.0]: https://github.com/HFooladi/chembl-pdb-linker/releases/tag/v0.1.0
