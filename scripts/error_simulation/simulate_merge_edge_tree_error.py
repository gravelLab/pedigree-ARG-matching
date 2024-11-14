import math
from itertools import combinations

from graph.coalescent_tree import CoalescentTree
from scripts.utility import *


def run_independent_simulations(tree: CoalescentTree):
    for simulation_counter in range(simulation_number):
        simulation_tree_name = f"{simulation_counter}"
        simulation_error_vertices = list(possible_error_vertices)
        for i in range(error_number):
            error_vertex = random.choice(simulation_error_vertices)
            simulation_error_vertices.remove(error_vertex)
            error_vertex_parent = tree.get_vertex_parent(error_vertex)
            coalescent_tree.merge_edge(parent=error_vertex_parent, child=error_vertex, recalculate_levels=False)
        tree.initialize_vertex_to_level_map()
        tree.save_to_file(simulation_tree_name)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)


def run_all_simulations(tree: CoalescentTree):
    for simulation_counter, errors_combination in enumerate(combinations(possible_error_vertices, error_number)):
        simulation_tree_name = f"{simulation_counter}"
        for error_vertex in errors_combination:
            error_vertex_parent = tree.get_vertex_parent(error_vertex)
            tree.merge_edge(parent=error_vertex_parent, child=error_vertex, recalculate_levels=False)
        tree.remove_isolated_vertices(recalculate_levels=True)
        tree.save_to_file(simulation_tree_name)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)


results_directory_paths = get_non_existing_directory_path("Specify the path where the results are to be stored:")
coalescent_tree_path = get_file_path("Specify the path to the coalescent tree:")
coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
possible_error_vertices = [x for x in coalescent_tree.children_map if x in coalescent_tree.parents_map]
max_error_number = len(possible_error_vertices)
error_number = get_natural_number_input_in_bounds(input_request=f"Specify the error number to be simulated"
                                                                f" (max value is {max_error_number}):",
                                                  lower_bound=1,
                                                  upper_bound=max_error_number)

error_mode_input = ("Specify the error-simulation mode:\n"
                    "1) Perform given number of independent simulations\n"
                    "2) Generate all the possible errors\n")

error_mode = get_natural_number_input_in_bounds(error_mode_input, 1, 2)
possible_trees_number = math.comb(max_error_number, error_number)
simulation_number = possible_trees_number
os.makedirs(results_directory_paths)
os.chdir(results_directory_paths)

if error_mode == 1:
    custom_simulation_number_input_request = (f"Specify the number of simulations. Notice that there are "
                                              f"{possible_trees_number} possible trees with {error_number}"
                                              f" errors\n"
                                              f"Any number higher than that will be treated by simulating "
                                              f"all the possible errors (the 2 option)")
    simulation_number = get_natural_number_input(custom_simulation_number_input_request)
    if simulation_number > possible_trees_number:
        error_mode = 2
        run_all_simulations(coalescent_tree)
    else:
        run_independent_simulations(coalescent_tree)
elif error_mode == 2:
    run_all_simulations(coalescent_tree)
