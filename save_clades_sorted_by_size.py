import os
from collections import defaultdict

from genealogical_graph import CoalescentTree, GenealogicalGraph
from utility import *

tree_filename = "clade"
pedigree_filename = "pedigree.pedigree"
proband_number_limit = 20
proband_number_examples_limit = 100
log_file_name = "log.txt"


def save_clades_by_proband_sizes(pedigree_filepath, pedigree: GenealogicalGraph):
    filepath = get_file_path("Specify the path to the coalescent tree:\n")
    result_filepath = get_non_existing_directory_name("Specify the path to save the clades to:\n")
    os.makedirs(result_filepath, exist_ok=False)
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath)
    # coalescent_tree.remove_unary_nodes()
    clades = coalescent_tree.get_connected_components()
    proband_number_to_clades = defaultdict(list)
    for clade in clades:
        # if len(clade) > 2 * proband_number_limit:
        #     continue
        proband_number = len(coalescent_tree.probands.intersection(clade))
        if proband_number > proband_number_limit:
            continue
        proband_number_to_clades[proband_number].append(clade)
    os.chdir(result_filepath)
    log_file = open(log_file_name, "w")
    log_file.write(f"The pedigree path: {pedigree_filepath}\n")
    log_file.write(f"The coalescent tree path {filepath}\n")
    log_file.close()
    for proband_number, clades in proband_number_to_clades.items():
        dir_name = str(proband_number)
        os.mkdir(dir_name)
        os.chdir(dir_name)
        for index, clade in enumerate(clades):
            if index > proband_number_examples_limit:
                break
            sub_dir_name = str(index)
            os.mkdir(sub_dir_name)
            os.chdir(sub_dir_name)
            coalescent_tree.get_subtree_from_vertices(clade).save_to_file(tree_filename)
            pedigree.save_ascending_genealogy_to_file(filename=pedigree_filename, probands=clade)
            os.chdir("..")
        os.chdir("..")


pedigree_filepath = get_file_path("Specify the path to the pedigree file:\n")
print("Processing the graph")
pedigree = GenealogicalGraph.get_diploid_graph_from_file(filename=pedigree_filepath)
save_clades_by_proband_sizes(pedigree_filepath, pedigree)
