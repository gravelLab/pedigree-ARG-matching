import pandas as pd
import numpy as np
import os

from scripts.utility.basic_utility import get_filepath, get_natural_number_input, get_non_existing_path


def bin_data(df, x_column, y_columns, bin_size):
    # Create bins based on the proband_number (x_column)
    bins = np.arange(df[x_column].min(), df[x_column].max() + bin_size, bin_size)
    df_binned = pd.cut(df[x_column], bins=bins, include_lowest=True, right=False)

    # Sum the values for each bin
    binned_df = df.groupby(df_binned, observed=False).agg({col: 'sum' for col in y_columns}).reset_index()

    # Extract the lowest bin boundary (left boundary) as the proband_number
    binned_df[x_column] = binned_df[x_column].apply(lambda x: x.left)

    # Reorder columns so x_column is first
    binned_df = binned_df[[x_column] + y_columns]
    return binned_df


def main():
    input_filepath = get_filepath("Specify the path to the data file:\n")
    # Derive output directory from input file
    input_dir = os.path.dirname(os.path.abspath(input_filepath))
    os.chdir(input_dir)
    output_filename = get_non_existing_path("Specify the output CSV"
                                            " filename (without extension): ").strip()
    if not output_filename.lower().endswith(".csv"):
        output_filename += ".csv"
    output_filepath = os.path.join(input_dir, output_filename)
    bin_size = get_natural_number_input("Specify the bin size:\n")
    df = pd.read_csv(input_filepath)
    x_column = 'proband_number'
    if x_column not in df.columns:
        raise ValueError(f"Error: The CSV file must contain a '{x_column}' column.")
    # Automatically detect numeric y_columns (excluding x_column)
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    y_columns = [col for col in numeric_columns if col != x_column]
    if not y_columns:
        raise ValueError("No numeric columns found to bin other than the x-column.")
    print(f"Detected y-columns for binning: {', '.join(y_columns)}")
    binned_df = bin_data(df, x_column, y_columns, bin_size)
    # Save binned data to output CSV
    binned_df.to_csv(output_filepath, index=False)
    print(f"Binned data saved to {output_filepath}")


if __name__ == "__main__":
    main()
