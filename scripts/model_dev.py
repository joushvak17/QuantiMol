import pandas as pd
from chembl_webresource_client.new_client import new_client


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

    return None


if __name__ == "__main__":
    main()
