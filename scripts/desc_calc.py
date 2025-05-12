import pandas as pd

from typing import List, Optional
from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    FindMolChiralCenters,
    GraphDescriptors,
    rdMolDescriptors as rdmd,
)


def compute_weiner_index(mol: Chem.Mol) -> Optional[float]:
    """Computes the Wiener index for a given molecule.

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        Optional[float]: The Wiener index of the molecule, or None if the molecule is invalid.
    """
    if mol is None:
        return None
    num_atoms = mol.GetNumAtoms()
    if num_atoms == 0:
        return None
    amat = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            res += amat[i][j]
    return res


def compute_descriptors(smiles_list: List[str]) -> pd.DataFrame:
    """Calculates various molecular descriptors for a list of SMILES strings using RDKit.

    Args:
        smiles_list (List[str]): A list of SMILES strings.

    Returns:
        pd.DataFrame: DataFrame containing the calculated molecular descriptors.
    """
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    descriptor_data = {
        "Molecular Weight": [Descriptors.ExactMolWt(mol) for mol in mol_list],
        "Number of Rotatable Bonds": [
            Descriptors.NumRotatableBonds(mol) for mol in mol_list
        ],
        "Number of Atoms": [mol.GetNumAtoms() for mol in mol_list],
        "Number of Bonds": [mol.GetNumBonds() for mol in mol_list],
        "Count of Chiral Centers": [
            len(FindMolChiralCenters(mol, includeUnassigned=True)) for mol in mol_list
        ],
        "Number of Rings": [rdmd.CalcNumRings(mol) for mol in mol_list],
        "Number of Aromatic Rings": [
            rdmd.CalcNumAromaticRings(mol) for mol in mol_list
        ],
        "Number of Hydrogen Bond Donors": [rdmd.CalcNumHBD(mol) for mol in mol_list],
        "Number of Hydrogen Bond Acceptors": [rdmd.CalcNumHBA(mol) for mol in mol_list],
        "Balaban J Index": [GraphDescriptors.BalabanJ(mol) for mol in mol_list],
        "Wiener Index": [compute_weiner_index(mol) for mol in mol_list],
        "LogP": [rdmd.CalcCrippenDescriptors(mol)[0] for mol in mol_list],
        "TPSA": [rdmd.CalcTPSA(mol) for mol in mol_list],
    }

    return pd.DataFrame(descriptor_data)
