from alignment.alignment_result import FailedClimbingCladeAlignmentResults
from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.configuration import logs_default_directory_name
from alignment.graph_matcher import GraphMatcher, get_initial_simulation_mapping_for_mode
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.alignment.run_alignment import save_alignment_result_to_files, get_store_and_save_vertex_alignment_callback
from scripts.utility.alignment_utility import read_mapping_from_file, read_initial_assignments_from_driver_file
from scripts.utility.basic_utility import *

alignment_general_results_filename = "results.txt"


def get_initial_mapping():
    initial_mapping_source_prompt = ("Specify how you want to get the initial mapping:\n"
                                     "1) Read from a text file\n"
                                     "2) Read from a driver file\n"
                                     "3) Get from the initial coalescent tree\n")
    initial_mapping_source_option = get_number_input_in_bounds(input_request=initial_mapping_source_prompt,
                                                               lower_bound=1,
                                                               upper_bound=3)
    match initial_mapping_source_option:
        case 1:
            initial_mapping_file = get_filepath("Specify the path to the initial mapping file:")
            return read_mapping_from_file(initial_mapping_file)
        case 2:
            driver_filepath = get_filepath("Specify the path to the driver file:")
            return read_initial_assignments_from_driver_file(driver_filepath)
        case 3:
            initial_coalescent_tree_filepath = get_filepath("Specify the path to the coalescent tree:")
            coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=initial_coalescent_tree_filepath)
            return get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)


def get_trees_sorted_by_index(trees_paths: str | Path):
    dir_path = Path(trees_paths)
    files = []
    for f in dir_path.iterdir():
        if not f.is_file():
            continue
        try:
            num = int(f.stem) if f.suffix else int(f.name)
            files.append((num, f.resolve()))
        except ValueError:
            continue
    return [path for _, path in sorted(files, key=lambda x: x[0])]


def run_alignments_with_all_trees(pedigree: PotentialMrcaProcessedGraph, initial_mapping: dict,
                                  result_directory: Path):
    alignment_number_with_solutions = []
    coalescent_trees_folder_files = [
        os.path.abspath(f) for f in os.listdir('.')
    ]
    for coalescent_tree_path in coalescent_trees_folder_files:
        filename = os.path.basename(coalescent_tree_path)
        # Ensure it's a non-metadata file
        if not os.path.isfile(coalescent_tree_path) or filename.startswith("."):
            continue
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
        results_directory = result_directory / filename
        alignment_general_results, save_results_callback = get_store_and_save_vertex_alignment_callback(
            coalescent_tree=coalescent_tree,
            directory_path=results_directory
        )
        graph_matcher = GraphMatcher(processed_graph=pedigree,
                                     coalescent_tree=coalescent_tree,
                                     initial_mapping=initial_mapping,
                                     logs_path=logs_default_directory_name,
                                     result_callback_function=save_results_callback)
        graph_matcher.find_alignments()
        save_alignment_result_to_files(graph_matcher=graph_matcher,
                                       alignment_result=alignment_general_results,
                                       directory_path=results_directory
                                       )
        clade_alignment_results = alignment_general_results.get_unique_clade_results()
        if isinstance(clade_alignment_results, FailedClimbingCladeAlignmentResults):
            continue
        if clade_alignment_results.alignments:
            alignment_number_with_solutions.append(filename)
    general_results_filepath = result_directory / alignment_general_results_filename
    alignment_general_results_file = open(general_results_filepath, 'w')
    alignment_general_results_file.write(f"Solutions are present in {alignment_number_with_solutions}\n")


def run_interactive_session():
    pedigree_path = get_filepath("Specify the absolute path to the pedigree:")
    coalescent_trees_folder = get_directory_path("Specify the path to the coalescent trees directory:")
    initial_mapping = get_initial_mapping()
    os.chdir(coalescent_trees_folder)
    result_directory = Path(get_non_existing_path("Specify the name of the results directory:"))
    os.mkdir(result_directory)
    pedigree_probands = {y for x in initial_mapping.values() for y in x}
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path,
                                                                         separation_symbol="\t",
                                                                         missing_parent_notation=["0"],
                                                                         preprocess_graph=True,
                                                                         probands=pedigree_probands,
                                                                         skip_first_line=True
                                                                         )
    run_alignments_with_all_trees(
        pedigree=pedigree,
        initial_mapping=initial_mapping,
        result_directory=result_directory
    )


if __name__ == "__main__":
    run_interactive_session()
