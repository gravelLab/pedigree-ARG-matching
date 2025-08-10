from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.PloidPedigree import PloidPedigree

from alignment.configuration import ProbandInitialAssignmentsMode
from alignment.graph_matcher import get_pedigree_simulation_probands_for_alignment_mode
from scripts.utility.basic_utility import *


def run_interactive_session():
    filepath = get_filepath("Specify the path to the coalescent tree. It should consist of one clade for more "
                            "meaningful results:\n")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    clade_number = len(coalescent_tree.get_connected_components())
    if clade_number != 1:
        raise ValueError(f"The coalescent tree must have 1 clade, found {clade_number} instead")
    tree_probands = coalescent_tree.get_sink_vertices()
    pedigree_probands = get_pedigree_simulation_probands_for_alignment_mode(
                                                                vertices=tree_probands,
                                                                alignment_mode=ProbandInitialAssignmentsMode.INDIVIDUAL
    )
    pedigree_filepath = get_filepath("Specify the path to the pedigree file:\n")
    genealogical_graph = PloidPedigree.get_ploid_pedigree_from_file(filepath=pedigree_filepath,
                                                                    probands=pedigree_probands)
    saving_format = get_natural_number_input_in_bounds("Specify the output format:\n"
                                                       "1) Tree-pedigree directories\n"
                                                       "2) Trees grouped by proband size\n", 1, 2)
    simulation_dir_name = get_non_existing_path("Specify the name for the simulation directory:\n")
    os.makedirs(simulation_dir_name)
    os.chdir(simulation_dir_name)
    probands = list(coalescent_tree.get_sink_vertices())
    tests_per_step = 100
    values_for_simulation = list(range(10, 100))

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
                coalescent_tree.save_ascending_graph_to_file(filepath="clade", vertices=random_probands)
                genealogical_graph.save_ascending_genealogy_as_diploid(filepath=f"{i}.pedigree",
                                                                       vertices=unphased_probands)
                os.chdir("..")
            else:
                coalescent_tree.save_ascending_graph_to_file(filepath=f"{i}", vertices=random_probands)
        os.chdir("..")


if __name__ == '__main__':
    run_interactive_session()
