from graph.coalescent_tree import CoalescentTree, SimpleGraph, GenealogicalGraph
from scripts.utility import *

filepath = get_file_path("Specify the path to the coalescent tree. It should consist of one clade for more "
                         "meaningful results:\n")
coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath)
pedigree_filepath = get_file_path("Specify the path to the pedigree file:\n")
simple_graph = SimpleGraph.get_diploid_graph_from_file(filepath=pedigree_filepath)
coalescent_tree_individuals = {x // 2 for x in coalescent_tree.probands}
unphased_probands = ([2 * x for x in coalescent_tree_individuals] +
                     [2 * x + 1 for x in coalescent_tree_individuals])
genealogical_graph = GenealogicalGraph(pedigree=simple_graph, probands=unphased_probands)
# print("Processing the graph")
# genealogical_graph = GenealogicalGraph(pedigree=Graph.get_pedigree_from_file(filename=pedigree_filepath),
#                                        probands=coalescent_tree.probands)
# print("The pedigree has been processed")
simulation_dir_name = get_non_existing_directory_path("Specify the name for the simulation directory:\n")
os.makedirs(simulation_dir_name)
os.chdir(simulation_dir_name)
proband_number = len(coalescent_tree.probands)
probands = list(coalescent_tree.probands)
reduction_step = 300
tests_per_step = 14
values_for_simulation = [value for value in range(7638, 100, -400) if value > 100]
values_for_simulation.extend([value for value in range(100, 9, -10) if value > 9])

for probands_left in values_for_simulation:
    os.makedirs(f"{probands_left}")
    os.chdir(f"{probands_left}")
    for i in range(tests_per_step):
        print(f"Running {i + 1} simulation")
        random_probands = random.sample(probands, probands_left)
        os.makedirs(f"{i}")
        os.chdir(f"{i}")
        coalescent_tree.save_ascending_genealogy_to_file(filename="clade", probands=random_probands)
        probands_corresponding_individuals = {x // 2 for x in random_probands}
        unphased_probands = ([2 * x for x in probands_corresponding_individuals] +
                             [2 * x + 1 for x in probands_corresponding_individuals])
        genealogical_graph.save_ascending_genealogy_to_file(filename=f"{i}.pedigree", probands=unphased_probands)
        os.chdir("..")
    os.chdir("..")
