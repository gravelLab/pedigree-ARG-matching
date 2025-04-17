from alignment.configuration import InitialMatchingMode
from alignment.graph_matcher import get_pedigree_simulation_probands_for_alignment_mode
from graph.coalescent_tree import CoalescentTree, SimpleGraph, GenealogicalGraph
from scripts.utility import *


def run_interactive_session():
    filepath = get_filepath("Specify the path to the coalescent tree. It should consist of one clade for more "
                            "meaningful results:\n")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    pedigree_filepath = get_filepath("Specify the path to the pedigree file:\n")
    simple_graph = SimpleGraph.get_diploid_graph_from_file(filepath=pedigree_filepath)
    pedigree_probands = get_pedigree_simulation_probands_for_alignment_mode(
                                                                coalescent_tree=coalescent_tree,
                                                                alignment_mode=InitialMatchingMode.INDIVIDUAL
    )
    saving_format = get_natural_number_input_in_bounds("Specify the output format:\n"
                                                       "1) Tree-pedigree directories\n"
                                                       "2) Trees grouped by proband size\n", 1, 2)
    genealogical_graph = GenealogicalGraph(pedigree=simple_graph, probands=pedigree_probands)
    simulation_dir_name = get_non_existing_path("Specify the name for the simulation directory:\n")
    os.makedirs(simulation_dir_name)
    os.chdir(simulation_dir_name)
    probands = list(coalescent_tree.probands)
    tests_per_step = 100
    values_for_simulation = list(range(40, 80))

    for probands_left in values_for_simulation:
        os.makedirs(f"{probands_left}")
        os.chdir(f"{probands_left}")
        for i in range(tests_per_step):
            print(f"Running {i + 1} simulation")
            random_probands = random.sample(probands, probands_left)
            probands_corresponding_individuals = {x // 2 for x in random_probands}
            unphased_probands = ([2 * x for x in probands_corresponding_individuals] +
                                 [2 * x + 1 for x in probands_corresponding_individuals])
            if saving_format == 1:
                os.makedirs(f"{i}")
                os.chdir(f"{i}")
                coalescent_tree.save_ascending_genealogy_to_file(filepath="clade", probands=random_probands)
                genealogical_graph.save_ascending_genealogy_to_file(filepath=f"{i}.pedigree", probands=unphased_probands)
                os.chdir("..")
            else:
                coalescent_tree.save_ascending_genealogy_to_file(filepath=f"{i}", probands=random_probands)
        os.chdir("..")


if __name__ == '__main__':
    run_interactive_session()
