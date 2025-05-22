import msprime.ancestry
from tskit import TreeSequence, Tree

from graph.coalescent_tree import CoalescentTree
from scripts.utility.basic_utility import *

sequence_length = 4


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
    tskit_tree.draw_svg(path=f"{coalescent_tree_file}.svg", size=(1000, 1000))


def run_interactive_session():
    # Asking the user for the path to the pedigree
    pedigree_path = get_filepath("Specify the path to the pedigree:")
    pedigree_file = open(pedigree_path, 'r')
    # Asking for the input parameters
    output_directory_path = input(
        "Specify the path to the result directory (where the simulated ARG will be saved):")
    if not os.path.exists(output_directory_path):
        os.makedirs(output_directory_path)

    pedigree = msprime.parse_pedigree(text_file=pedigree_file, sequence_length=sequence_length)

    individuals_dict = dict()
    for id, ind in enumerate(pedigree.individuals):
        mapped_value = int(ind.metadata["file_id"])
        individuals_dict[2 * id] = 2 * mapped_value
        individuals_dict[2 * id + 1] = 2 * mapped_value + 1

    pedigree = pedigree.tree_sequence()

    simulated_arg = msprime.ancestry.sim_ancestry(initial_state=pedigree,
                                                  model="fixed_pedigree",
                                                  sequence_length=sequence_length)
    simulated_arg: TreeSequence
    counter = 0
    for tree in simulated_arg.trees():
        coalescent_tree_filename = f"{output_directory_path}/coalescent_tree_{counter}"
        save_coalescent_tree(tree, individuals_dict, coalescent_tree_filename)
        counter += 1


if __name__ == '__main__':
    run_interactive_session()
