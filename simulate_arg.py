from utility import *

import msprime.ancestry
from tskit import TreeSequence, Tree
from genealogical_graph import CoalescentTree


def save_coalescent_tree(tskit_tree: Tree, individual_id_dict: dict, filename: str):
    coalescent_tree_file = open(filename, 'w')
    coalescent_tree = CoalescentTree.get_coalescent_tree(tskit_tree)
    coalescent_tree: CoalescentTree

    for level in coalescent_tree.levels:
        for vertex in level:
            if vertex in coalescent_tree.parents_map:
                [parent] = coalescent_tree.parents_map[vertex]
                parent = individual_id_dict[parent]
                vertex = individual_id_dict[vertex]
                coalescent_tree_file.write(f"{vertex} {parent}\n")
    coalescent_tree_file.close()


# Asking the user for the path to the pedigree
pedigree_path = get_file_path("Specify the path to the pedigree:")
pedigree_file = open(pedigree_path, 'r')
# Asking for the input parameters
output_directory_path = get_file_path(
    "Specify the path to the result directory (where the simulated ARG will be saved):")

pedigree = msprime.parse_pedigree(text_file=pedigree_file, sequence_length=1)

individuals_dict = dict()
for id, ind in enumerate(pedigree.individuals):
    mapped_value = int(ind.metadata["file_id"])
    individuals_dict[2 * id] = 2 * mapped_value
    individuals_dict[2 * id + 1] = 2 * mapped_value + 1

pedigree = pedigree.tree_sequence()

simulated_arg = msprime.ancestry.sim_ancestry(initial_state=pedigree,
                                              model="fixed_pedigree",
                                              sequence_length=1)
simulated_arg: TreeSequence
counter = 0
for tree in simulated_arg.trees():
    coalescent_tree_filename = f"{output_directory_path}/coalescent_tree_{counter}"
    save_coalescent_tree(tree, individuals_dict, coalescent_tree_filename)
    counter += 1
