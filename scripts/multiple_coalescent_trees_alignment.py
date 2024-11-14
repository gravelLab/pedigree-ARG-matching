import os

from alignment.graph_matcher import MatcherLogger, GraphMatcher
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.coalescent_tree import CoalescentTree
from scripts.run_alignment import save_alignment_result_to_files
from scripts.utility import get_file_path, get_directory_path, read_mapping_from_file, get_non_existing_directory_path

alignment_general_results_filename = "results.txt"

pedigree_path = get_file_path("Specify the absolute path to the pedigree:")
coalescent_trees_folder = get_directory_path("Specify the path to the coalescent trees directory:")
initial_mapping_file = get_file_path("Specify the path to the initial mapping:")
os.chdir(coalescent_trees_folder)
result_directory = get_non_existing_directory_path("Specify the name of the results directory:")
os.mkdir(result_directory)
initial_mapping = read_mapping_from_file(initial_mapping_file)
pedigree_probands = {y for x in initial_mapping.values() for y in x}
pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path,
                                                                     separation_symbol="\t",
                                                                     missing_parent_notation=["0"],
                                                                     initialize_levels=False,
                                                                     initialize_ancestor_maps=False)
pedigree.reduce_to_ascending_genealogy(pedigree_probands, recalculate_levels=True)
pedigree.initialize_potential_mrca_map()

alignment_number_with_solutions = []
current_index = 0
coalescent_trees_folder_files = list(os.listdir('.'))
for coalescent_tree_filename in coalescent_trees_folder_files:
    if os.path.isfile(coalescent_tree_filename):  # Ensure it's a file
        coalescent_tree_path = os.path.abspath(coalescent_tree_filename)
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(coalescent_tree_path)
        os.chdir(result_directory)
        logger = MatcherLogger()
        graph_matcher = GraphMatcher(processed_graph=pedigree,
                                     coalescent_tree=coalescent_tree,
                                     initial_mapping=initial_mapping,
                                     logger=logger)

        dir_name = f"{current_index}"
        os.mkdir(dir_name)
        os.chdir(dir_name)
        alignment_results = graph_matcher.find_mapping()
        save_alignment_result_to_files(alignment_results, coalescent_tree, pedigree)
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
