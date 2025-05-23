import os
import joblib
import pandas as pd
import lightgbm as lgb

from loguru import logger
from typing import List
from chembl_webresource_client.new_client import new_client
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from molfeat.trans.fp import FPVecTransformer

from descriptors import compute_descriptors as descriptors


def main() -> None:
    """Main function to fetch, process, and train data from ChEMBL.

    Returns:
        None: Returns None after processing the data.
    """

    # Define constants
    IMPORTANT_COLUMNS: List = [
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_value",
    ]
    UPPER_THRESHOLD: int = 1000

    # Create an activity client
    activity_client = new_client.activity

    # Prompt the user for a search string
    search_string: str = input("Enter a search string: ")  # e.g., "Cancer"
    logger.info(f"Searching for: {search_string}")

    # Search for the given string
    activity_query = activity_client.search(search_string)

    # Convert to a DataFrame and print the number of records
    activity_df = pd.DataFrame.from_dict(activity_query)
    if activity_df.empty:
        logger.error("No records found.")
        return None
    else:
        logger.info(f"Number of records found: {len(activity_df)}")

    # Filter the DataFrame
    activity_df = (
        activity_df[IMPORTANT_COLUMNS]
        .drop_duplicates(["canonical_smiles"])
        .dropna(subset=["canonical_smiles", "standard_value"])
        .astype({"standard_value": float})
        .query("standard_value <= 100 or standard_value >= @UPPER_THRESHOLD")
    )

    # Append the classification column and print the number of records
    activity_df["class"] = activity_df["standard_value"].apply(
        lambda x: "inactive" if x >= UPPER_THRESHOLD else "active"
    )
    logger.info(f"Number of records after filtering: {len(activity_df)}")

    # Calculate descriptors
    descriptor_df = descriptors(activity_df["canonical_smiles"].tolist())

    # Concatenate the DataFrames
    df_combined = pd.concat([activity_df, descriptor_df], axis=1)

    # Define the desired column order
    insert_after_column = 1
    original_cols = activity_df.columns.tolist()
    new_order = (
        original_cols[: insert_after_column + 1]  # Before insertion point
        + descriptor_df.columns.tolist()  # Descriptor columns
        + original_cols[insert_after_column + 1 :]  # After insertion point
    )

    df = df_combined[new_order]

    # Drop the "standard_value" column and start with the third column
    # df = df.drop(df[-2], axis=1, inplace=True)
    df = df.drop(columns=["standard_value"])
    df = df.iloc[:, 2:]

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
    logger.info(f"Train size: {X_train.shape[0]}")
    logger.info(f"Test size: {X_test.shape[0]}")
    logger.info(f"Number of features: {X_train.shape[1]}")
    logger.info(f"Number of classes: {len(label_encoder.classes_)}")

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

    # Save the model to the models directory
    model_path = os.path.join("src", "models", "model.pkl")
    transformer.to_state_yaml_file(model_path.replace(".pkl", ".yaml"))
    joblib.dump(model, model_path)
    logger.info(f"Model saved to {model_path}")

    return None


if __name__ == "__main__":
    main()
