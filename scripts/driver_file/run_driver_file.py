import argparse

from alignment.driver_file import ProcessedDriverFile
from alignment.graph_matcher import *
from scripts.alignment.run_alignment import (run_alignments,
                                             get_store_and_save_vertex_alignment_callback,
                                             save_alignment_result_to_files)
from scripts.utility.basic_utility import *


def run_driver_file(driver_filepath):
    processed_driver_file = ProcessedDriverFile.process_driver_file(filepath=driver_filepath)
    processed_driver_file.preprocess_graphs_for_alignment()
    alignment_general_results, save_results_callback = get_store_and_save_vertex_alignment_callback(
        coalescent_tree=processed_driver_file.coalescent_tree,
        directory_path=processed_driver_file.output_path
    )
    graph_matcher = GraphMatcher(coalescent_tree=processed_driver_file.coalescent_tree,
                                 processed_graph=processed_driver_file.pedigree,
                                 initial_mapping=processed_driver_file.initial_assignments,
                                 logs_path=processed_driver_file.output_path,
                                 alignment_vertex_mode=processed_driver_file.alignment_vertex_mode,
                                 alignment_edge_mode=processed_driver_file.alignment_edge_mode,
                                 result_callback_function=save_results_callback
                                 )
    run_alignments(graph_matcher=graph_matcher)
    save_alignment_result_to_files(
        alignment_result=alignment_general_results,
        coalescent_tree=processed_driver_file.coalescent_tree,
        pedigree=processed_driver_file.pedigree,
        directory_path=processed_driver_file.output_path
    )


def main():
    parser = argparse.ArgumentParser(description="Run graph alignment.")
    parser.add_argument("--file", "-f", type=str, help="Path to the driver file")
    args = parser.parse_args()
    if args.file:
        driver_filepath = args.file
    else:
        driver_filepath = get_filepath("Specify the path to the driver file:")

    run_driver_file(driver_filepath=driver_filepath)


if __name__ == "__main__":
    main()
