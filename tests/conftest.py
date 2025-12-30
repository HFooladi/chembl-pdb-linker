"""Shared pytest fixtures for ChEMBL-PDB Linker tests."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest
from chembl_pdb_linker.config import Config


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def config(temp_dir):
    """Create a test configuration with temp directories."""
    config = Config(base_dir=temp_dir)
    return config


@pytest.fixture
def sample_chembl_activities():
    """Sample ChEMBL activity data."""
    return pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL25", "CHEMBL59", "CHEMBL1201", "CHEMBL941"],
            "canonical_smiles": [
                "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
                "CC(C)Cc1ccc(C(C)C(=O)O)cc1",  # Ibuprofen
                "Cc1ccc(C(=O)c2ccccc2)cc1",  # Benzophenone
                "CC(=O)Nc1ccc(O)cc1",  # Paracetamol
            ],
            "standard_inchi_key": [
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "HEFNNWSXXWATRW-UHFFFAOYSA-N",
                "SDDLEVPIDBLVHC-UHFFFAOYSA-N",
                "RZVAJINKPMORJF-UHFFFAOYSA-N",
            ],
            "uniprot_id": ["P23219", "P23219", "P35354", "P23219"],
            "target_chembl_id": ["CHEMBL221", "CHEMBL221", "CHEMBL210", "CHEMBL221"],
            "target_name": [
                "Cyclooxygenase-1",
                "Cyclooxygenase-1",
                "Cannabinoid receptor 1",
                "Cyclooxygenase-1",
            ],
            "standard_type": ["IC50", "IC50", "Ki", "IC50"],
            "standard_value": [100.0, 50.0, 1500.0, 200.0],
            "standard_units": ["nM", "nM", "nM", "nM"],
            "pchembl_value": [7.0, 7.3, 5.8, 6.7],
            "confidence_score": [9, 9, 9, 9],
            "assay_chembl_id": ["CHEMBL1000001", "CHEMBL1000002", "CHEMBL1000003", "CHEMBL1000004"],
        }
    )


@pytest.fixture
def sample_sifts_mapping():
    """Sample SIFTS PDB-UniProt mapping data."""
    return pd.DataFrame(
        {
            "pdb_id": ["1PRH", "1PRH", "3PGH", "4U0J", "5KIR"],
            "chain_id": ["A", "B", "A", "A", "A"],
            "uniprot_id": ["P23219", "P23219", "P35354", "P23219", "P23219"],
            "res_begin": [1, 1, 10, 1, 1],
            "res_end": [300, 300, 350, 300, 300],
            "pdb_begin": [1, 1, 10, 1, 1],
            "pdb_end": [300, 300, 350, 300, 300],
        }
    )


@pytest.fixture
def sample_ligand_data():
    """Sample PDB ligand data with InChIKeys."""
    return pd.DataFrame(
        {
            "pdb_id": ["1PRH", "3PGH", "4U0J"],
            "ligand_code": ["ASP", "IBU", "SAL"],
            "ligand_inchikey": [
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",  # Matches aspirin
                "HEFNNWSXXWATRW-UHFFFAOYSA-N",  # Matches ibuprofen
                "YGSDEFSMJLZEOE-UHFFFAOYSA-N",  # Different compound
            ],
            "ligand_name": ["Aspirin", "Ibuprofen", "Salicylic acid"],
        }
    )


@pytest.fixture
def sample_linked_data(sample_chembl_activities, sample_sifts_mapping):
    """Sample linked ChEMBL-PDB data."""
    # Join on uniprot_id
    linked = sample_chembl_activities.merge(
        sample_sifts_mapping[["pdb_id", "chain_id", "uniprot_id"]].drop_duplicates(),
        on="uniprot_id",
        how="inner",
    )
    # Group PDB IDs per activity
    grouped = (
        linked.groupby(list(sample_chembl_activities.columns))
        .agg({"pdb_id": list, "chain_id": list})
        .reset_index()
    )
    grouped = grouped.rename(columns={"pdb_id": "pdb_ids", "chain_id": "chain_ids"})
    return grouped


@pytest.fixture
def sample_activities_mixed_units():
    """Sample activities with different unit types for conversion testing."""
    return pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4", "CHEMBL5"],
            "standard_value": [100.0, 1.0, 0.001, 10.0, 500.0],
            "standard_units": ["nM", "uM", "mM", "pM", "nM"],
            "standard_type": ["IC50", "IC50", "IC50", "Ki", "Kd"],
            "confidence_score": [9, 9, 8, 9, 7],
        }
    )


@pytest.fixture
def sample_pdb_metadata():
    """Sample PDB entry metadata."""
    return pd.DataFrame(
        {
            "pdb_id": ["1PRH", "3PGH", "4U0J", "5KIR", "2ABC"],
            "title": [
                "Structure of COX-1",
                "Structure of CB1",
                "Complex of COX-1 with inhibitor",
                "High resolution COX-1",
                "NMR structure of protein",
            ],
            "resolution": [2.5, 2.8, 1.9, 1.5, None],  # None for NMR
            "experimental_method": [
                "X-ray diffraction",
                "X-ray diffraction",
                "X-ray diffraction",
                "X-ray diffraction",
                "NMR",
            ],
            "release_date": ["2000-01-01", "2005-06-15", "2014-08-20", "2016-03-10", "2010-12-01"],
        }
    )


@pytest.fixture
def real_output_path():
    """Path to the real pipeline output for integration tests."""
    path = Path(__file__).parent.parent / "data" / "curated" / "bioactivity_pdb_linked.parquet"
    if not path.exists():
        pytest.skip("Real pipeline output not available. Run the pipeline first.")

    # Verify this is the real output, not test sample data
    df = pd.read_parquet(path)
    if len(df) < 1000:
        pytest.skip("Output file contains sample data, not real pipeline output.")

    return path
