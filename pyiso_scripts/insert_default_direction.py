import pandas as pd
import argparse

def insert_weights(weights_file, features_file, output_file):
    # Read the weights file as a DataFrame
    weights_df = pd.read_csv(weights_file, sep='\t', comment='#', header=None, skip_blank_lines=True)

    # Extract feature names and normalized weights from the weights DataFrame
    feature_names = weights_df.iloc[0].tolist()  # The first row contains feature names
    normalized_weights = weights_df.iloc[1].tolist()  # The second row contains normalized weights

    # Create a mapping of feature names to their weights
    weights_map = dict(zip(feature_names, normalized_weights))

    # Read the features file as a DataFrame
    features_df = pd.read_csv(features_file, sep='\t')

    # Construct the insert row dynamically based on feature matching
    insert_row = ["DefaultDirection"]
    for col in features_df.columns[1:]:
        if col in weights_map:
            insert_row.append(weights_map[col])  # Insert corresponding weight
        else:
            insert_row.append("-")  # Fill with "-" for unmatched columns

    # Insert the new row into the DataFrame
    insert_df = pd.DataFrame([insert_row], columns=features_df.columns)
    updated_df = pd.concat([insert_df, features_df]).reset_index(drop=True)
    updated_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert weights into features file.")
    parser.add_argument("--weights_pin", required=True, help="Path to the weights.pin file generated from a normal Percolator run.")
    parser.add_argument("--features_pin", required=True, help="Path to the input make-pin.pin file.")
    parser.add_argument("--output_pin", required=True, help="Path to save the inserted.features.pin file with normalized final output weights inserted as the default direction.")
    args = parser.parse_args()

    try:
        insert_weights(args.weights_pin, args.features_pin, args.output_pin)
        print(f"[INFO]Default direction successfully inserted into {args.output_pin}.")
    except Exception as e:
        print(f"[Error]{e}")
