from graph.coalescent_tree import CoalescentTree, GenealogicalGraph
from scripts.utility.basic_utility import *


def save_largest_clade_and_get_probands(input_filepath: str | Path, output_filepath: str | Path):
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=input_filepath)
    largest_clade = coalescent_tree.get_subtree_from_vertices(coalescent_tree.get_largest_clade_by_probands())
    largest_clade.save_to_file(output_filepath)
    return set(largest_clade.probands)


def run_interactive_session():
    tree_input_filepath = get_filepath("Specify the path to the coalescent tree:\n")
    tree_output_filepath = get_filepath("Specify the path to the clade resulting file:\n")
    pedigree_input_filepath = get_filepath("Specify the path to the pedigree file:\n")
    pedigree_output_filepath = input("Specify the path to the resulting file for the ascending genealogy:\n")
    print("Processing the tree")
    probands = save_largest_clade_and_get_probands(input_filepath=tree_input_filepath,
                                                   output_filepath=tree_output_filepath)
    print("Processing the graph")
    genealogical_graph = GenealogicalGraph.get_diploid_graph_from_file(filepath=pedigree_input_filepath)
    genealogical_graph.save_ascending_genealogy_to_file(filepath=pedigree_output_filepath, probands=probands)
