import os
import warnings
from functools import partial

from alignment.alignment_result import AlignmentResult, FailedClimbingAlignmentResult, FullAlignmentResult
from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.configuration import AlignmentVertexMode, AlignmentEdgeMode, logs_default_directory_name
from alignment.graph_matcher import GraphMatcher, get_initial_simulation_mapping_for_mode
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.utility.basic_utility import get_filepath, get_directory_path, get_non_existing_path

results_directory_name = "results"
results_filename = "proband_number_to_alignments_number.csv"


def collect_vertex_alignment_number(alignment_number_counter: int, alignment_result: AlignmentResult):
    match alignment_result:
        case FailedClimbingAlignmentResult(clade_root=clade_root, failed_vertex=failed_vertex):
            raise ValueError(f"Failed climbing for the clade with root {clade_root} at {failed_vertex}")
        case FullAlignmentResult(vertex_alignment=vertex_alignment):
            assert vertex_alignment
            alignment_number_counter += 1


def run_interactive_session():
    pedigree_filepath = get_filepath("Specify the path to the pedigree:")
    clades_filepath = get_directory_path("Specify the path to the directory containing the subclades grouped by "
                                         "proband number:")
    # Create the directories to save the results
    results_directory_path = os.path.join(clades_filepath, results_directory_name)
    os.makedirs(results_directory_path, exist_ok=True)
    current_directory = os.getcwd()
    os.chdir(results_directory_path)
    simulation_name = get_non_existing_path("Specify the simulation name:")
    simulation_results_path = os.path.join(results_directory_path, simulation_name)
    simulation_logs_path = os.path.join(results_directory_path, simulation_name, logs_default_directory_name)
    results_filepath = os.path.join(simulation_results_path, results_filename)
    os.mkdir(simulation_results_path)
    os.mkdir(simulation_logs_path)
    os.chdir(current_directory)
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_filepath)
    # Run the alignments
    with open(results_filepath, 'a') as results_file:
        for proband_number_directory in os.listdir(clades_filepath):
            try:
                proband_number = int(proband_number_directory)
            except ValueError:
                warnings.warn("The directory's name is not a valid proband number, skipping")
                continue
            proband_number_directory_path = os.path.join(clades_filepath, proband_number_directory)
            for tree_file in os.listdir(proband_number_directory_path):
                tree_filepath = os.path.join(proband_number_directory_path, tree_file)
                coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_filepath)
                coalescent_tree.remove_unary_nodes()
                initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
                alignment_counter = 0
                callback_function = partial(collect_vertex_alignment_number, alignment_counter)
                log_directory_path = os.path.join(simulation_logs_path, proband_number_directory, tree_file)
                os.makedirs(log_directory_path, exist_ok=True)
                matcher = GraphMatcher(coalescent_tree=coalescent_tree,
                                       processed_graph=pedigree,
                                       logs_path=log_directory_path,
                                       initial_mapping=initial_mapping,
                                       alignment_vertex_mode=AlignmentVertexMode.ALL_ALIGNMENTS,
                                       alignment_edge_mode=AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT,
                                       result_callback_function=callback_function)
                matcher.find_alignments()
                # Save the results
                results_file.write(f"{proband_number},{alignment_counter}")


if __name__ == '__main__':
    run_interactive_session()
