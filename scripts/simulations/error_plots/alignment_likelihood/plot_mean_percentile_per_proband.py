import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from scripts.utility.basic_utility import get_filepath


def interactive_session():
    # Load the CSV with averages
    input_filepath = get_filepath("Specify the path to the csv file with the mean data: ")
    input_filepath = Path(input_filepath)
    parent_directory = input_filepath.parent
    os.chdir(parent_directory)
    output_filename = input("Enter the filename of the output file (press Enter to use "
                            "the same filename as the input file):")
    if not output_filename:
        output_filename = input_filepath.stem
    if not output_filename.endswith('.svg'):
        output_filename = f"{output_filename}.svg"
    output_filepath = parent_directory / output_filename
    df = pd.read_csv(input_filepath)
    # Plot
    plt.figure(figsize=(10, 6))
    plt.bar(df["proband_number"].astype(str), df["identity_percentile"], color='skyblue')
    # plt.plot(df["proband_number"], df["identity_percentile"], marker='o', linestyle='-')
    plt.xlabel("Proband Number")
    plt.ylabel("Average Identity Percentile")
    plt.title("Average Identity Percentile per Proband")
    # plt.grid(True)
    # plt.tight_layout()
    # Save or show
    plt.savefig(output_filepath, dpi=300)
    plt.show()


if __name__ == "__main__":
    interactive_session()
