import pandas as pd
import numpy as np

from scripts.utility import get_filepath, get_non_existing_path, get_natural_number_input


def bin_data(df, x_column, y_columns, bin_size):
    # Create bins based on the proband_number (x_column)
    bins = np.arange(df[x_column].min(), df[x_column].max() + bin_size, bin_size)
    df_binned = pd.cut(df[x_column], bins=bins, include_lowest=True, right=False)
    # Sum the values for each bin
    binned_df = df.groupby(df_binned, observed=False).agg({col: 'sum' for col in y_columns}).reset_index()
    # Extract the lowest bin boundary (left boundary) as the proband_number
    binned_df[x_column] = binned_df[x_column].apply(lambda x: x.left)
    # Reorder the columns so proband_number is the first column
    binned_df = binned_df[[x_column] + y_columns]
    return binned_df


def main():
    input_filepath = get_filepath("Specify the path to the data file:\n")
    output_filepath = get_non_existing_path("Specify the output CSV file path:\n")
    bin_size = get_natural_number_input("Specify the bin size:\n")
    # Read input CSV file
    df = pd.read_csv(input_filepath)
    # Define x-axis and values to bin
    x_column = 'proband_number'
    y_columns = ['no_solutions', 'individual_and_spouses', 'individual_and_non_spouse',
                 'no_individual_spouse', 'neither_individual_nor_spouse']
    # Perform binning
    binned_df = bin_data(df, x_column, y_columns, bin_size)
    # Save binned data to a new CSV file
    binned_df.to_csv(output_filepath, index=False)
    print(f"Binned data saved to {output_filepath}")


if __name__ == "__main__":
    main()
