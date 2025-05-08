import os
from pathlib import Path

import pandas as pd

from scripts.utility import get_filepath, get_non_existing_path


def group_and_save_the_average(input_filepath: str | Path, output_filepath: str | Path):
    df = pd.read_csv(input_filepath)
    # Group by 'proband_number' and compute the mean
    avg_per_proband = (
        df.groupby("proband_number")["identity_percentile"]
        .mean()
        .reset_index()
    )
    try:
        avg_per_proband["proband_number"] = avg_per_proband["proband_number"].astype(int)
        avg_per_proband = avg_per_proband.sort_values("proband_number")
    except ValueError:
        # If conversion fails (some values aren't integers), skip sorting
        pass
    avg_per_proband.to_csv(output_filepath, index=False)


def interactive_session():
    input_filepath = get_filepath("Enter the filepath of the input file: ")
    input_filepath = Path(input_filepath)
    parent_directory = input_filepath.parent
    os.chdir(parent_directory)
    output_filename = get_non_existing_path("Enter the filename of the output file: ")
    output_filepath = parent_directory / output_filename
    group_and_save_the_average(input_filepath=input_filepath, output_filepath=output_filepath)


if __name__ == "__main__":
    interactive_session()
