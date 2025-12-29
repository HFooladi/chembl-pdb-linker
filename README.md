# ChEMBL-PDB Linker

Link ChEMBL bioactivity data with PDB structural information to create a curated dataset of compounds with both activity measurements and 3D co-crystal structures.

## Features

- **Protein-level linking**: Connect ChEMBL targets to PDB structures via UniProt IDs
- **Ligand-level linking**: Match ChEMBL compounds to PDB ligands via InChIKey
- **Configurable filters**: Activity types (IC50, Ki, Kd, etc.), confidence scores, resolution
- **Hybrid data fetching**: Bulk downloads for main data, API for cross-references
- **Parquet output**: Efficient, compressed output format

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

## Usage

### Command Line

```bash
# Run complete pipeline
chembl-pdb-linker run

# Or run steps separately:
chembl-pdb-linker download          # Download ChEMBL and PDB data
chembl-pdb-linker link              # Link datasets via UniProt/InChIKey
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

## Configuration

Edit `config/default.yaml` to customize:

```yaml
chembl:
  confidence_score: 9          # Filter for single-protein targets
  activity_types:              # Activity types to include
    - IC50
    - Ki
    - Kd
  standard_units:
    - nM

pdb:
  max_resolution: 3.5          # Maximum structure resolution (Å)

linking:
  use_inchikey: true           # Match ligands by InChIKey
  inchikey_connectivity_only: false  # Use full InChIKey or just first 14 chars
```

## Output Schema

The final Parquet file contains:

| Column | Description |
|--------|-------------|
| `chembl_id` | ChEMBL compound ID |
| `smiles` | Canonical SMILES |
| `inchikey` | Standard InChIKey |
| `uniprot_id` | Target UniProt accession |
| `target_name` | Target protein name |
| `activity_type` | IC50/Ki/Kd/etc. |
| `activity_value_nm` | Activity value in nM |
| `pchembl` | pChEMBL value |
| `pdb_id` | PDB structure ID |
| `pdb_ligand_code` | 3-letter HET code |
| `resolution` | Structure resolution (Å) |
| `rcsb_url` | Download URL for structure |

## Data Sources

- **ChEMBL**: [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/)
- **PDB**: [https://www.rcsb.org/](https://www.rcsb.org/)
- **SIFTS**: [https://www.ebi.ac.uk/pdbe/docs/sifts/](https://www.ebi.ac.uk/pdbe/docs/sifts/)

## License

MIT
