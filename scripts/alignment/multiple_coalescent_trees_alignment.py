from alignment.configuration import logs_default_directory_name
from alignment.graph_matcher import GraphMatcher, get_initial_simulation_mapping_for_mode
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.coalescent_tree import CoalescentTree
from scripts.alignment.run_alignment import save_alignment_result_to_files
from scripts.utility.basic_utility import *
from scripts.utility.alignment_utility import read_mapping_from_file

alignment_general_results_filename = "results.txt"


def run_interactive_session():
    pedigree_path = get_filepath("Specify the absolute path to the pedigree:")
    coalescent_trees_folder = get_directory_path("Specify the path to the coalescent trees directory:")
    initial_mapping_source_prompt = ("Specify how you want to get the initial mapping:\n"
                                     "1) Read from a file\n"
                                     "2) Get from the initial coalescent tree\n")
    initial_mapping_source_option = get_natural_number_input_in_bounds(input_request=initial_mapping_source_prompt,
                                                                       lower_bound=1,
                                                                       upper_bound=2)
    if initial_mapping_source_option == 1:
        initial_mapping_file = get_filepath("Specify the path to the initial mapping:")
        initial_mapping = read_mapping_from_file(initial_mapping_file)
    else:
        initial_coalescent_tree_filepath = get_filepath("Specify the path to the coalescent tree:")
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=initial_coalescent_tree_filepath)
        initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
        del coalescent_tree
    os.chdir(coalescent_trees_folder)
    result_directory = get_non_existing_path("Specify the name of the results directory:")
    os.mkdir(result_directory)
    pedigree_probands = {y for x in initial_mapping.values() for y in x}
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path,
                                                                         separation_symbol="\t",
                                                                         missing_parent_notation=["0"],
                                                                         preprocess_graph=False
                                                                         )
    pedigree.reduce_to_ascending_genealogy(pedigree_probands, recalculate_levels=True)
    pedigree.initialize_potential_mrca_map()

    alignment_number_with_solutions = []
    current_index = 0
    coalescent_trees_folder_files = list(os.listdir('.'))
    for coalescent_tree_filename in coalescent_trees_folder_files:
        # Ensure it's a file
        if os.path.isfile(coalescent_tree_filename):
            coalescent_tree_path = os.path.abspath(coalescent_tree_filename)
            coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
            os.chdir(result_directory)
            graph_matcher = GraphMatcher(processed_graph=pedigree,
                                         coalescent_tree=coalescent_tree,
                                         initial_mapping=initial_mapping,
                                         logs_path=logs_default_directory_name)

            dir_name = f"{current_index}"
            os.mkdir(dir_name)
            os.chdir(dir_name)
            alignment_results = graph_matcher.find_mapping()
            save_alignment_result_to_files(pedigree=pedigree,
                                           coalescent_tree=coalescent_tree,
                                           alignment_result=alignment_results
                                           )
            clade_root = next(iter(alignment_results))
            alignments = alignment_results[clade_root]
            if alignments:
                alignment_number_with_solutions.append(current_index)
            current_index += 1
            os.chdir("..")
            os.chdir("..")

    os.chdir(result_directory)
    alignment_general_results_file = open(alignment_general_results_filename, 'w')
    alignment_general_results_file.write(f"Solutions are present in {alignment_number_with_solutions}\n")


if __name__ == "__main__":
    run_interactive_session()
