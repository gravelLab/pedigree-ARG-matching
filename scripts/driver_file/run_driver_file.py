import argparse

from alignment.driver_file import ProcessedDriverFile
from alignment.graph_matcher import *
from scripts.alignment.run_alignment import run_alignments_and_save_results
from scripts.utility.basic_utility import *


def main():
    parser = argparse.ArgumentParser(description="Run graph alignment.")
    parser.add_argument("--file", "-f", type=str, help="Path to the driver file")
    args = parser.parse_args()
    if args.file:
        driver_filepath = args.file
    else:
        driver_filepath = get_filepath("Specify the path to the driver file:")

    processed_driver_file = ProcessedDriverFile.process_driver_file(filepath=driver_filepath)
    processed_driver_file.preprocess_graphs_for_alignment()
    graph_matcher = GraphMatcher(coalescent_tree=processed_driver_file.coalescent_tree,
                                 processed_graph=processed_driver_file.pedigree,
                                 initial_mapping=processed_driver_file.initial_assignments,
                                 logs_path=processed_driver_file.output_path,
                                 matching_mode=processed_driver_file.alignment_mode)
    run_alignments_and_save_results(graph_matcher=graph_matcher, output_path=processed_driver_file.output_path)


if __name__ == "__main__":
    main()
