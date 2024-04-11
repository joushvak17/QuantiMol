from multiprocessing import Pool

import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import (rdMolDescriptors as
                        rdmd,
                        GraphDescriptors,
                        Descriptors,
                        FindMolChiralCenters)


def compute_descriptors(smile):
    mol = Chem.MolFromSmiles(smile)

    amat = Chem.GetDistanceMatrix(mol)
    res = 0
    num_atoms = mol.GetNumAtoms()
    for j in range(num_atoms):
        for k in range(j + 1, num_atoms):
            res += amat[j][k]
    wiener_res = res

    descriptor_data = {
        'Molecular Weight': Descriptors.ExactMolWt(mol),
        'Number of Rotatable Bonds': rdmd.CalcNumRotatableBonds(mol),
        'Number of Atoms': mol.GetNumAtoms(),
        'Number of Bonds': mol.GetNumBonds(),
        'Count of Chiral Centers': len(FindMolChiralCenters(mol, includeUnassigned=True)),
        'Number of Rings': rdmd.CalcNumRings(mol),
        'Number of Aromatic Rings': rdmd.CalcNumAromaticRings(mol),
        'Number of Hydrogen Bond Donors': rdmd.CalcNumHBD(mol),
        'Number of Hydrogen Bond Acceptors': rdmd.CalcNumHBA(mol),
        'Balaban J Index': GraphDescriptors.BalabanJ(mol),
        'Wiener Index': wiener_res,
        'LogP': rdmd.CalcCrippenDescriptors(mol)[0],
        'TPSA': rdmd.CalcTPSA(mol),
    }

    return descriptor_data


def descriptors_calculation(smiles_list, compute_function, num_processes):
    os.environ['OMP_NUM_THREADS'] = '1'
    with Pool(processes=num_processes) as pool:
        results = pool.map(compute_function, smiles_list)

    descriptor_values = pd.DataFrame(results)

    return descriptor_values.round(2)
