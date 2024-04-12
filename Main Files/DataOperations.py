import pandas as pd
from rdkit import Chem
from rdkit.Chem import (rdMolDescriptors as
                        rdmd,
                        GraphDescriptors,
                        Descriptors,
                        FindMolChiralCenters)


def descriptors_calculation(smiles_list):
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    wiener_res = []
    amat_list = [Chem.GetDistanceMatrix(mol) for mol in mol_list]
    for i, mol in enumerate(mol_list):
        res = 0
        amat = amat_list[i]
        num_atoms = mol.GetNumAtoms()
        for j in range(num_atoms):
            for k in range(j + 1, num_atoms):
                res += amat[j][k]
        wiener_res.append(res)

    descriptor_data = {
        'Molecular Weight': [Descriptors.ExactMolWt(mol) for mol in mol_list],
        'Number of Rotatable Bonds': [rdmd.CalcNumRotatableBonds(mol) for mol in mol_list],
        'Number of Atoms': [mol.GetNumAtoms() for mol in mol_list],
        'Number of Bonds': [mol.GetNumBonds() for mol in mol_list],
        'Count of Chiral Centers': [len(FindMolChiralCenters(mol, includeUnassigned=True)) for mol in
                                    mol_list],
        'Number of Rings': [rdmd.CalcNumRings(mol) for mol in mol_list],
        'Number of Aromatic Rings': [rdmd.CalcNumAromaticRings(mol) for mol in mol_list],
        'Number of Hydrogen Bond Donors': [rdmd.CalcNumHBD(mol) for mol in mol_list],
        'Number of Hydrogen Bond Acceptors': [rdmd.CalcNumHBA(mol) for mol in mol_list],
        'Balaban J Index': [GraphDescriptors.BalabanJ(mol) for mol in mol_list],
        'Wiener Index': wiener_res,
        'LogP': [rdmd.CalcCrippenDescriptors(mol)[0] for mol in mol_list],
        'TPSA': [rdmd.CalcTPSA(mol) for mol in mol_list],
    }

    descriptor_values = pd.DataFrame(descriptor_data)

    return descriptor_values.round(2)
