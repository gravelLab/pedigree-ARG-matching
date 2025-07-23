import os
from pathlib import Path

import pandas as pd

from scripts.utility.basic_utility import get_filepath, get_non_existing_path, get_natural_number_with_lower_bound


def bin_simulation_results(input_csv, output_csv, bin_size):
    df = pd.read_csv(input_csv)
    # Ignore 'proband_number' and 'row_count'
    classification_columns = df.columns[2:]
    binned_results = []
    temp_rows = []
    accumulated_count = 0
    for _, row in df.iterrows():
        temp_rows.append(row)
        accumulated_count += int(row['row_count'])
        if accumulated_count >= bin_size:
            # Aggregate the bin
            temp_df = pd.DataFrame(temp_rows)
            # Compute weighted mean for classification columns
            weighted_averages = (temp_df[classification_columns].multiply(temp_df['row_count'], axis=0).sum()
                                 / temp_df['row_count'].sum())
            binned_results.append([temp_df['proband_number'].iloc[0], accumulated_count] + list(weighted_averages))
            temp_rows = []
            accumulated_count = 0
    binned_df = pd.DataFrame(binned_results, columns=['proband_number', 'row_count'] + list(classification_columns))
    binned_df[classification_columns] = binned_df[classification_columns].round(3)
    binned_df['proband_number'] = binned_df['proband_number'].astype(int)
    binned_df['row_count'] = binned_df['row_count'].astype(int)
    binned_df.to_csv(output_csv, index=False)
    print(f"Binned results saved to {output_csv}")


if __name__ == "__main__":
    input_filepath = Path(get_filepath("Specify the path to the input file:"))
    os.chdir(input_filepath.parent)
    output_filepath = get_non_existing_path("Specify the output filename (without extension):")
    output_filepath = f"{output_filepath}.csv"
    input_bin_size = get_natural_number_with_lower_bound("Specify the bin size:", lower_bound=2)
    bin_simulation_results(input_filepath, output_filepath, input_bin_size)
