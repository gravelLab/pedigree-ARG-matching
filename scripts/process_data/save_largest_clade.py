from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.PloidPedigree import PloidPedigree

from scripts.utility.basic_utility import *


def save_largest_clade_and_get_probands(input_filepath: str | Path, output_filepath: str | Path):
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=input_filepath)
    coalescent_tree.reduce_to_subgraph(coalescent_tree.get_largest_clade_by_probands())
    coalescent_tree.save_to_file(output_filepath)
    return set(coalescent_tree.get_sink_vertices())


def run_interactive_session():
    tree_input_filepath = get_filepath("Specify the path to the coalescent tree:\n")
    tree_output_filepath = get_non_existing_path("Specify the path to the clade resulting file:\n")
    pedigree_input_filepath = get_filepath("Specify the path to the pedigree file:\n")
    pedigree_output_filepath = get_non_existing_path("Specify the path to the resulting"
                                                     " file for the ascending pedigree:\n")
    print("Processing the tree")
    probands = save_largest_clade_and_get_probands(input_filepath=tree_input_filepath,
                                                   output_filepath=tree_output_filepath)
    print(f"The largest clade has {len(probands)} probands")
    print("Processing the graph")
    genealogical_graph = PloidPedigree.get_ploid_pedigree_from_file(filepath=pedigree_input_filepath)
    genealogical_graph.save_ascending_genealogy_as_diploid(filepath=pedigree_output_filepath, vertices=probands)


if __name__ == '__main__':
    run_interactive_session()
