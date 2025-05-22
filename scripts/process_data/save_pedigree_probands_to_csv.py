import csv

from graph.genealogical_graph import GenealogicalGraph
from scripts.utility.basic_utility import get_filepath, get_non_existing_path


def run_interactive_session():
    pedigree_filepath = get_filepath("Specify the path to the pedigree file:")
    result_filepath = get_non_existing_path("Specify the path to the output csv file (without extension):")
    result_filepath = f"{result_filepath}.csv"
    pedigree = GenealogicalGraph.get_diploid_graph_from_file(filepath=pedigree_filepath)
    probands = pedigree.get_probands()
    with open(result_filepath, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write all probands in one row
        writer.writerow(probands)


if __name__ == '__main__':
    run_interactive_session()
