import os
from pathlib import Path

import pandas as pd

from scripts.utility import get_filepath, get_non_existing_path


def run_interactive_session():
    filepath = Path(get_filepath("Specify the path to the data file:"))
    os.chdir(filepath.parent)
    output_filepath = get_non_existing_path("Specify the output filename (without extension):")
    output_filepath = f"{output_filepath}.csv"
    df = pd.read_csv(filepath)
    classification_columns = df.columns[1:]

    # Compute row-wise sums and normalize
    row_sums = df[classification_columns].sum(axis=1)
    df[classification_columns] = df[classification_columns].div(row_sums, axis=0)

    # Group by proband_number and calculate mean for each class, also count the rows in each group
    grouped_df = df.groupby('proband_number').agg(
        {col: 'mean' for col in classification_columns},  # Calculate mean for all classification columns
    ).reset_index()

    # Count the number of rows in each group
    count_df = df.groupby('proband_number').size().reset_index(name='row_count')

    # Merge the count_df with grouped_df on 'proband_number'
    grouped_df = pd.merge(grouped_df, count_df, on='proband_number', how='left')

    # Rearrange the columns to make 'row_count' the second column
    columns_order = ['proband_number', 'row_count'] + [col for col in classification_columns]
    grouped_df = grouped_df[columns_order]

    # Round the classification values to 3 decimal places
    grouped_df[classification_columns] = grouped_df[classification_columns].round(3)
    grouped_df.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    run_interactive_session()
