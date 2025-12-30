# ChEMBL-PDB Linker

Link ChEMBL bioactivity data with PDB structural information to create a curated dataset of protein-ligand pairs with both activity measurements and 3D co-crystal structures. Similar to PDBbind but derived from ChEMBL and PDB directly.

## Features

- **Validated protein-ligand pairs**: Ensures PDB structure contains BOTH the target protein AND the ligand
- **Protein-level linking**: Connect ChEMBL targets to PDB structures via UniProt IDs
- **Ligand-level linking**: Match ChEMBL compounds to PDB ligands via InChIKey
- **RCSB Search API integration**: Efficient bulk queries for ligand-to-structure mappings
- **Configurable filters**: Activity types (IC50, Ki, Kd), confidence scores, resolution
- **Fully reproducible**: End-to-end pipeline from raw data to curated output
- **Parquet output**: Efficient, compressed output format

## Expected Output

Running the full pipeline produces:
- **~98,500** validated protein-ligand pairs with bioactivity data
- **~9,000** unique compounds
- **~14,700** unique PDB structures  
- **~1,300** unique target proteins

## Installation

```bash
# Clone the repository
git clone https://github.com/HFooladi/chembl-pdb-linker.git
cd chembl-pdb-linker

# Option 1: Install with uv (recommended)
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .

# Option 2: Install with pip
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Dependencies

- Python >= 3.9
- pandas, pyarrow
- rdkit
- httpx, requests
- typer, pyyaml, tqdm

## Quick Start

```bash
# Run complete pipeline (downloads ~5GB, takes ~30-60 min)
chembl-pdb-linker run
```

The output is saved to `data/curated/bioactivity_pdb_linked.parquet`.

## Usage

### Command Line

```bash
# Run complete pipeline
chembl-pdb-linker run

# Or run steps separately:
chembl-pdb-linker download          # Download ChEMBL, SIFTS, and PDB ligand data
chembl-pdb-linker link              # Link datasets with protein-ligand validation
chembl-pdb-linker extract           # Generate final curated dataset

# Show statistics
chembl-pdb-linker stats

# Use custom config
chembl-pdb-linker run --config my_config.yaml
```

### Python API

```python
from chembl_pdb_linker import Pipeline, Config

# Load config and create pipeline
config = Config.default()
pipeline = Pipeline(config)

# Run complete pipeline
output_path = pipeline.run()

# Or run steps separately
pipeline.download(chembl_version="34")
linked_df = pipeline.link()
output_path = pipeline.extract()

# Get statistics
stats = pipeline.get_statistics()
```

## How It Works

### Pipeline Overview

```
ChEMBL Database ──┬─→ Extract activities (confidence=9) ──┐
                  │                                        │
SIFTS Mapping ────┼─→ UniProt ↔ PDB mapping ──────────────┼──→ Protein-level linking
                  │                                        │
PDB Ligands ──────┼─→ Ligand code ↔ InChIKey mapping ─────┤
                  │                                        │
RCSB Search API ──┴─→ Ligand code → PDB structures ───────┴──→ Validated pairs
```

### Key Algorithm: Protein-Ligand Pair Validation

The critical step is ensuring each linked pair is **validated**:

1. **Protein matching**: ChEMBL target → UniProt ID → PDB structures (via SIFTS)
2. **Ligand matching**: ChEMBL compound → InChIKey → PDB ligand codes
3. **Validation**: For each potential match, verify the PDB structure contains BOTH the protein AND the ligand


## Configuration

Edit `config/default.yaml` to customize:

```yaml
chembl:
  confidence_score: 9          # Filter for single-protein targets
  activity_types:              # Activity types to include
    - IC50
    - Ki
    - Kd
    - EC50
  standard_units:
    - nM

pdb:
  max_resolution: 3.5          # Maximum structure resolution (Å)
  rcsb_search_api: "https://search.rcsb.org/rcsbsearch/v2/query"

linking:
  use_inchikey: true           # Match ligands by InChIKey
  inchikey_connectivity_only: false  # Use full InChIKey or just first 14 chars
```

## Output Schema

The final Parquet file (`data/curated/bioactivity_pdb_linked.parquet`) contains:

| Column | Description |
|--------|-------------|
| `chembl_id` | ChEMBL compound ID |
| `smiles` | Canonical SMILES |
| `inchikey` | Standard InChIKey |
| `uniprot_id` | Target UniProt accession |
| `target_name` | Target protein name |
| `activity_type` | IC50/Ki/Kd/EC50 |
| `activity_value` | Activity value |
| `activity_unit` | Unit (typically nM) |
| `pchembl` | pChEMBL value (-log10(M)) |
| `pdb_id` | PDB structure ID |
| `pdb_ligand_code` | 3-letter HET code |
| `resolution` | Structure resolution (Å) |
| `rcsb_url` | Download URL for structure |

## Data Sources

- **ChEMBL**: [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/) - Bioactivity data
- **PDB/RCSB**: [https://www.rcsb.org/](https://www.rcsb.org/) - 3D structures
- **SIFTS**: [https://www.ebi.ac.uk/pdbe/docs/sifts/](https://www.ebi.ac.uk/pdbe/docs/sifts/) - UniProt-PDB mappings
- **Ligand Expo**: [http://ligand-expo.rcsb.org/](http://ligand-expo.rcsb.org/) - Ligand InChIKey data

## Comparison to PDBbind

| Aspect | PDBbind | ChEMBL-PDB Linker |
|--------|---------|-------------------|
| Size | ~20K complexes | ~98K pairs |
| Affinity source | Literature curation | ChEMBL database |
| Structure source | PDB | PDB |
| Update frequency | Annual | On-demand |
| Reproducibility | Manual | Fully automated |

## License

MIT
