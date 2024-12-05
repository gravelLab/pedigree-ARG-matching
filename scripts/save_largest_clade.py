from graph.coalescent_tree import CoalescentTree, GenealogicalGraph
from scripts.utility import *


def save_largest_clade_and_get_probands():
    filepath = get_filepath("Specify the path to the coalescent tree:\n")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    largest_clade = coalescent_tree.get_subtree_from_vertices(coalescent_tree.get_largest_clade_by_probands())
    while True:
        try:
            result_filepath = input("Specify the path to the resulting file:\n")
            largest_clade.save_to_file(result_filepath)
            break
        except OSError:
            print("Specify a different file path")
    return set(largest_clade.probands)


probands = save_largest_clade_and_get_probands()
pedigree_filepath = get_filepath("Specify the path to the pedigree file:\n")
print("Processing the graph")
# genealogical_graph = GenealogicalGraph(pedigree=SimpleGraph.get_pedigree_from_file(filename=pedigree_filepath),
#                                        probands=probands)
genealogical_graph = GenealogicalGraph.get_graph_from_file(filepath=pedigree_filepath, ploidy=2)
result_filepath = input("Specify the path to the resulting file for the ascending genealogy:\n")
print("Saving the ascending genealogy")
# genealogical_graph.save_to_file(filename=result_filepath)
