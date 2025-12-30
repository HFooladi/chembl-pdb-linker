# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```bash
# Install (with uv)
uv venv && source .venv/bin/activate && uv pip install -e .

# Install dev dependencies
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Lint
ruff check src/

# Format code
ruff format src/

# Type check
pyright

# Run tests
pytest
```

## CLI Usage

```bash
# Run complete pipeline
chembl-pdb-linker run

# Run individual steps
chembl-pdb-linker download          # Download ChEMBL SQLite + SIFTS mapping + PDB ligands
chembl-pdb-linker link              # Link via UniProt/InChIKey with validation
chembl-pdb-linker extract           # Generate final parquet

# Show dataset statistics
chembl-pdb-linker stats

# Use custom config
chembl-pdb-linker run --config config/custom.yaml
```

## Architecture

This package links ChEMBL bioactivity data with PDB structural information to create a PDBbind-like dataset with validated protein-ligand pairs.

### Pipeline Flow (`src/chembl_pdb_linker/pipeline.py`)

**Download Phase:**
1. Download ChEMBL SQLite → Extract activities/compounds/targets
2. Download SIFTS mapping → Parse UniProt-PDB mappings
3. Load ligand InChIKey mapping (from `data/raw/pdb_ligand_inchikeys.parquet`)
4. Fetch PDB ligand-structure mapping via RCSB Search API

**Link Phase:**
5. Protein-level linking via UniProt IDs
6. Ligand-level linking via InChIKey with **protein-ligand pair validation**

**Extract Phase:**
7. Standardize bioactivity values
8. Enrich with structure metadata
9. Save final curated dataset

Data flows: `data/raw/` → `data/intermediate/` → `data/curated/`

### Core Components
1. **Download** → **Link** → **Extract**
2. `Pipeline` class orchestrates all components
3. Data flows through `data/raw/` → `data/intermediate/` → `data/curated/`

**Downloaders** (`src/chembl_pdb_linker/downloaders/`):
- `ChEMBLDownloader`: Downloads ChEMBL SQLite, extracts activities/compounds/targets via SQL queries
- `PDBDownloader`: Downloads SIFTS mapping (UniProt↔PDB), fetches PDB ligand-to-structure mappings via RCSB Search API

**Linkers** (`src/chembl_pdb_linker/linkers/`):
- `UniProtLinker`: Protein-level linking - matches ChEMBL targets to PDB structures via shared UniProt IDs
- `LigandLinker`: Compound-level linking with **validation** - ensures the PDB structure contains BOTH the target protein AND the ligand (not just matching InChIKeys independently)

**Extractors** (`src/chembl_pdb_linker/extractors/`):
- `BioactivityExtractor`: Standardizes activity units, computes pChEMBL values
- `StructureExtractor`: Enriches with PDB metadata, filters by resolution

### Key Implementation Details

**Protein-Ligand Pair Validation** (`linkers/ligand.py`):
- `_validate_protein_ligand_pair()`: Ensures `pdb_id` from ligand matching is in `pdb_ids` from protein matching
- Critical for avoiding false positives where ligand appears in a different PDB than the protein
- Typical filtering: ~99% of naive InChIKey matches are false positives

**RCSB Search API Integration** (`downloaders/pdb.py`):
- `fetch_ligand_structures_via_rcsb()`: Efficiently queries which PDB structures contain specific ligands
- Uses reverse lookup (query by ligand code) instead of per-PDB queries
- Required for large-scale datasets (82K+ PDB structures)

### Configuration (`src/chembl_pdb_linker/config.py`)
Dataclass-based config loaded from `config/default.yaml`. Key settings:
- `chembl.confidence_score`: Filter for single-protein targets (default: 9)
- `chembl.activity_types`: IC50, Ki, Kd, EC50
- `pdb.max_resolution`: Structure quality threshold (default: 3.5 Å)
- `pdb.rcsb_search_api`: RCSB Search API URL for ligand queries
- `linking.inchikey_connectivity_only`: Use first 14 chars of InChIKey for looser matching

### Data Files
- **Input**: ChEMBL SQLite, SIFTS CSV, `pdb_ligand_inchikeys.parquet` (ligand code → InChIKey)
- **Intermediate**: `chembl_activities.parquet`, `pdb_sifts_mapping.parquet`, `pdb_ligands.parquet`
- **Output**: `data/curated/bioactivity_pdb_linked.parquet`

### Expected Results
- ~98,500 validated protein-ligand pairs with bioactivity data
- ~9,000 unique compounds, ~14,700 PDB structures, ~1,300 proteins
- Similar to PDBbind but with ChEMBL bioactivity data
