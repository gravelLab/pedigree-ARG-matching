from collections import defaultdict

from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.PloidPedigree import PloidPedigree

from alignment.configuration import ProbandInitialAssignmentsMode
from alignment.graph_matcher import get_pedigree_simulation_probands_for_alignment_mode
from scripts.utility.basic_utility import *

tree_filename = "clade"
clade_ascending_pedigree_filename = "pedigree.pedigree"
proband_number_limit = 200
proband_number_tree_examples_limit = 100
info_filename = "info.txt"


def save_clades_by_proband_sizes(pedigree_filepath: str, coalescent_tree_filepath: str,
                                 result_filepath: str):
    print("Processing the graph")
    pedigree = PloidPedigree.get_ploid_pedigree_from_file(filepath=pedigree_filepath)
    print("Processing the coalescent tree")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_filepath)
    os.makedirs(result_filepath, exist_ok=False)
    clades = coalescent_tree.get_connected_components()
    proband_number_to_clades = defaultdict(list)
    clades_found = 0
    # Calculating the maximum number of examples that we can find. This can help us stop the search without processing
    # all the clades if we have already found everything we need
    # Specifically, for every proband number from [2, proband_number_limit], we will have at maximum
    # proband_number_tree_examples_limit clade-pedigree pairs
    clades_limit = (proband_number_limit - 1) * proband_number_tree_examples_limit
    sink_vertices = set(coalescent_tree.get_sink_vertices())
    for clade in clades:
        # Since we haven't removed the unary nodes, we can't estimate the maximum number of probands based on the
        # size of the clade
        proband_number = len(sink_vertices.intersection(clade))
        if proband_number < 2:
            continue
        if proband_number > proband_number_limit:
            continue
        if len(proband_number_to_clades[proband_number]) >= proband_number_tree_examples_limit:
            continue
        proband_number_to_clades[proband_number].append(clade)
        clades_found += 1
        if clades_found >= clades_limit:
            break
    os.chdir(result_filepath)
    log_file = open(info_filename, "w")
    log_file.write(f"The pedigree path: {pedigree_filepath}\n")
    log_file.write(f"The coalescent tree path {coalescent_tree_filepath}\n")
    log_file.close()
    for proband_number, clades in proband_number_to_clades.items():
        dir_name = str(proband_number)
        os.mkdir(dir_name)
        os.chdir(dir_name)
        for index, clade in enumerate(clades):
            sub_dir_name = str(index)
            os.mkdir(sub_dir_name)
            os.chdir(sub_dir_name)
            tree_copy = coalescent_tree.copy()
            tree_copy.reduce_to_subgraph(subgraph_vertices=clade)
            clade_probands = tree_copy.get_sink_vertices()
            tree_copy.save_to_file(filepath="clade")
            pedigree_probands = get_pedigree_simulation_probands_for_alignment_mode(
                vertices=clade_probands, alignment_mode=ProbandInitialAssignmentsMode.INDIVIDUAL
            )
            pedigree.save_ascending_genealogy_as_diploid(filepath=clade_ascending_pedigree_filename,
                                                         vertices=pedigree_probands)
            os.chdir("..")
        os.chdir("..")


def run_interactive_session():
    pedigree_filepath = get_filepath("Specify the path to the pedigree file:\n")
    coalescent_tree_filepath = get_filepath("Specify the path to the coalescent tree:\n")
    result_filepath = get_non_existing_path("Specify the path to save the clades to:\n")
    save_clades_by_proband_sizes(pedigree_filepath=pedigree_filepath, coalescent_tree_filepath=coalescent_tree_filepath,
                                 result_filepath=result_filepath)


if __name__ == '__main__':
    run_interactive_session()
