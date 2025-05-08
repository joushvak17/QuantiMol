import pandas as pd
from chembl_webresource_client.new_client import new_client

from rdkit import Chem
from rdkit.Chem import (
    rdMolDescriptors as rdmd,
    GraphDescriptors,
    Descriptors,
    FindMolChiralCenters,
)


def descriptors(smiles_list: list[str]) -> pd.DataFrame:
    """Calculates various molecular descriptors for a list of SMILES strings using RDKit.

    Args:
        smiles_list (list[str]): A list of SMILES strings representing the molecules.

    Returns:
        pd.DataFrame: DataFrame containing the calculated molecular descriptors.
    """
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
        "Wiener Index": wiener_res,
        "LogP": [rdmd.CalcCrippenDescriptors(mol)[0] for mol in mol_list],
        "TPSA": [rdmd.CalcTPSA(mol) for mol in mol_list],
    }

    descriptor_values = pd.DataFrame(descriptor_data)

    return descriptor_values


def main() -> None:
    """Main function to fetch, process, and train data from ChEMBL.

    Returns:
        None: Returns None after processing the data.
    """

    # Prompt user for a search string
    search_string: str = input("Enter a search string: ")
    print(f"Searching for: {search_string}")

    # Create a activity client and search for the given string
    activity_client = new_client.activity
    activity_query = activity_client.search(search_string)

    # Convert to a DataFrame and print the number of records
    activity_df = pd.DataFrame.from_records(activity_query)
    print(f"Number of records found: {len(activity_df)}")

    # Filter out NA values in the "standard_value" column
    activity_df = activity_df[activity_df.standard_value.notna()]

    # Filter out NA values in the "canonical_smiles" column
    activity_df = activity_df[activity_df.canonical_smiles.notna()]

    # Drop duplicates based on "canonical_smiles"
    activity_df = activity_df.drop_duplicates(["canonical_smiles"])

    # Set the important columns
    imp_col = ["molecule_chembl_id", "canonical_smiles", "standard_value"]
    activity_df = activity_df[imp_col]

    # Change the data type of "standard_value" to float
    activity_df["standard_value"] = activity_df["standard_value"].astype(float)

    # Filter out values between 100 and 1000
    activity_df = activity_df[~activity_df["standard_value"].between(100, 1000)]

    # Append the classification column
    activity_df["class"] = activity_df["standard_value"].apply(
        lambda x: "inactive" if x >= 1000 else "active"
    )

    # Calculate descriptors
    descriptor_df = descriptors(activity_df["canonical_smiles"].tolist())

    # Concatenate the descriptors with the activity data, placing the class column last
    insert_after_column: int = 1

    df_data_before = activity_df.iloc[:, : insert_after_column + 1]
    df_data_after = activity_df.iloc[:, insert_after_column + 1 :]
    df = pd.concat([df_data_before, descriptor_df, df_data_after], axis=1)

    return None


if __name__ == "__main__":
    main()
