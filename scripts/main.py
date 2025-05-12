import pandas as pd
import lightgbm as lgb
import joblib

from loguru import logger
from chembl_webresource_client.new_client import new_client
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from molfeat.trans.fp import FPVecTransformer


# TODO: Move this function to a separate file

# TODO: Add logging rather than print statements
logger.


def main() -> None:
    """Main function to fetch, process, and train data from ChEMBL.

    Returns:
        None: Returns None after processing the data.
    """

    # Prompt user for a search string
    search_string: str = input("Enter a search string: ")  # Cancer
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

    # Drop the "standard_value" column and start with the third column
    # df = df.drop(df[-2], axis=1, inplace=True)
    df = df.drop(columns=["standard_value"])
    df = df.iloc[:, 2:]

    # TODO: Move the model training to a separate file
    # Create the x and y data
    X = df.drop(columns=["class"])
    y = df["class"]

    # Create the transformer and transform the data
    smiles = activity_df["canonical_smiles"]
    transformer = FPVecTransformer(kind="desc2D", dtype=float)
    features = transformer(smiles)
    features_df = pd.DataFrame(features, index=smiles.index)
    X = pd.concat([X, features_df], axis=1)

    # Create the label encoder and transform the labels
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y)

    # Create the train and test data, splitting 80% for training and 20% for testing
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    print("Train size: " + str(X_train.shape[0]))
    print("Test size: " + str(X_test.shape[0]))
    print("Number of features: " + str(X_train.shape[1]))
    print("Number of classes: " + str(len(label_encoder.classes_)))

    # Create a LightGBM dataset
    train_data = lgb.Dataset(X_train, label=y_train)
    test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

    # Create the params
    params = {
        "task": "train",
        "objective": "binary",
        "boosting_type": "dart",
        "data_sample_strategy": "bagging",
        "tree_learner": "data",
        "metric": "binary_logloss",
        "num_leaves": 71,
        "learning_rate": 0.3,
        "feature_fraction": 0.9,
    }

    # Train the model
    model = lgb.train(
        params=params,
        train_set=train_data,
        num_boost_round=400,
        valid_sets=test_data,
    )

    # Make some predictions
    y_pred = model.predict(X_test, num_iteration=model.best_iteration)
    y_pred_binary = [1 if x > 0.5 else 0 for x in y_pred]
    accuracy = accuracy_score(y_test, y_pred_binary)
    print(f"Accuracy: {accuracy}")

    # Save the model
    # transformer.to_state_yaml_file()
    # joblib.dump(model, "model.pkl")

    return None


if __name__ == "__main__":
    main()
