import math
from itertools import combinations

from graph.coalescent_tree import CoalescentTree
from scripts.utility import *


def run_independent_simulations(tree: CoalescentTree, possible_error_vertices: [int], simulation_number: int,
                                merge_error_operations_per_simulation: int, coalescent_tree_path: str):
    for simulation_counter in range(simulation_number):
        simulation_tree_name = f"{simulation_counter}"
        simulation_error_vertices = list(possible_error_vertices)
        for i in range(merge_error_operations_per_simulation):
            error_vertex = random.choice(simulation_error_vertices)
            simulation_error_vertices.remove(error_vertex)
            error_vertex_parent = tree.get_vertex_parent(error_vertex)
            tree.merge_edge(parent=error_vertex_parent, child=error_vertex, recalculate_levels=False)
        tree.initialize_vertex_to_level_map()
        tree.save_to_file(simulation_tree_name)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)


def run_all_simulations(tree: CoalescentTree, possible_error_vertices: [int],
                        merge_error_operations_per_simulation: int, coalescent_tree_path: str):
    for simulation_counter, errors_combination in enumerate(combinations(possible_error_vertices,
                                                                         merge_error_operations_per_simulation)):
        simulation_tree_name = f"{simulation_counter}"
        for error_vertex in errors_combination:
            error_vertex_parent = tree.get_vertex_parent(error_vertex)
            tree.merge_edge(parent=error_vertex_parent, child=error_vertex, recalculate_levels=False)
        tree.remove_isolated_vertices(recalculate_levels=True)
        tree.save_to_file(simulation_tree_name)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)


def run_interactive_session_single_directory():
    results_directory_path = get_non_existing_path("Specify the path where the results are to be stored:")
    coalescent_tree_path = get_filepath("Specify the path to the coalescent tree:")
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
    possible_error_vertices = [x for x in coalescent_tree.children_map if x in coalescent_tree.parents_map]
    max_error_number = len(possible_error_vertices)
    merge_error_operations_per_simulation = get_natural_number_input_in_bounds(
        input_request=f"Specify the number of merge operations that"
                      f"are to be done per simulation"
                      f" (max value is {max_error_number}):",
        lower_bound=1,
        upper_bound=max_error_number)

    possible_trees_number = math.comb(max_error_number, merge_error_operations_per_simulation)
    simulation_mode_input = ("Specify the error-simulation mode:\n"
                             "1) Perform given number of simulations (with possible duplicates)\n"
                             f"2) Generate all the possible errors. There are {possible_trees_number} possible trees\n")

    simulation_mode = get_natural_number_input_in_bounds(simulation_mode_input, 1, 2)
    simulation_number = -1
    if simulation_mode == 1:
        custom_simulation_number_input_request = (f"Specify the number of simulations. Notice that there are "
                                                  f"{possible_trees_number} possible trees with "
                                                  f"{merge_error_operations_per_simulation} errors\n"
                                                  f"Any number higher than that will be treated by simulating "
                                                  f"all the possible errors (the 2 option)")
        simulation_number = get_natural_number_input(custom_simulation_number_input_request)
        if simulation_number > possible_trees_number:
            print("Generating all the possible modified trees")
            simulation_mode = 2

    os.makedirs(results_directory_path)
    os.chdir(results_directory_path)

    if simulation_mode == 1:
        run_independent_simulations(tree=coalescent_tree, coalescent_tree_path=coalescent_tree_path,
                                    possible_error_vertices=possible_error_vertices,
                                    merge_error_operations_per_simulation=merge_error_operations_per_simulation,
                                    simulation_number=simulation_number)
    elif simulation_mode == 2:
        run_all_simulations(tree=coalescent_tree, coalescent_tree_path=coalescent_tree_path,
                            merge_error_operations_per_simulation=merge_error_operations_per_simulation,
                            possible_error_vertices=possible_error_vertices)


def run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(
        parent_directory: Path,
        results_path: Path,
        merge_error_operations_per_simulation: int):
    for tree_pedigree_directory in os.listdir(parent_directory):
        tree_pedigree_directory_path = parent_directory / tree_pedigree_directory
        if not os.path.isdir(tree_pedigree_directory_path):
            continue
        paths = get_paths_from_tree_pedigree_directory(tree_pedigree_directory_path)
        if not paths:
            continue
        _, tree_path = paths
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path)
        possible_error_vertices = [x for x in coalescent_tree.children_map if x in coalescent_tree.parents_map]
        results_subpath = results_path / tree_pedigree_directory
        os.makedirs(results_subpath)
        os.chdir(results_subpath)
        run_all_simulations(tree=coalescent_tree, coalescent_tree_path=tree_path,
                            merge_error_operations_per_simulation=merge_error_operations_per_simulation,
                            possible_error_vertices=possible_error_vertices)
        os.chdir(parent_directory)


def run_interactive_session_multiple_tree_pedigree_subdirectories_all_trees():
    parent_directory = get_directory_path("Specify the path to the parent directory with tree-pedigree directories:")
    parent_directory = Path(parent_directory)
    results_directory_path = get_non_existing_path("Specify the path where the results are to be stored:")
    results_directory_path = Path(results_directory_path)
    merge_error_operations_per_simulation = get_natural_number_with_lower_bound(
        input_request=f"Specify the number of merge operations that will be used for all the simulations.\nIf the "
                      f"number of merge operations exceeds the possible number of errors, nothing will be done:",
        lower_bound=1
    )
    run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(
        parent_directory=parent_directory,
        results_path=results_directory_path,
        merge_error_operations_per_simulation=merge_error_operations_per_simulation
    )


def run_all_simulations_for_multiple_parent_directory_with_tree_pedigree_subdirectories(
        parent_directory: Path,
        results_path: Path,
        merge_error_operations_per_simulation: int):
    for tree_pedigrees_parent_directory in os.listdir(parent_directory):
        results_subpath = results_path / tree_pedigrees_parent_directory
        parent_subdirectory = parent_directory / tree_pedigrees_parent_directory
        if not os.path.isdir(parent_subdirectory):
            continue
        run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(
            parent_directory=parent_subdirectory,
            results_path=results_subpath,
            merge_error_operations_per_simulation=merge_error_operations_per_simulation
        )


def run_interactive_session_multiple_parent_directories_with_tree_pedigree_subdirectories_all_trees():
    parent_directory = get_directory_path("Specify the path to the parent directory with subdirectories containing"
                                          " tree-pedigree directories:")
    parent_directory = Path(parent_directory)
    results_directory_path = get_non_existing_path("Specify the path where the results are to be stored:")
    results_directory_path = Path(results_directory_path)
    merge_error_operations_per_simulation = get_natural_number_with_lower_bound(
        input_request=f"Specify the number of merge operations that will be used for all the simulations.\nIf the "
                      f"number of merge operations exceeds the possible number of errors, nothing will be done:",
        lower_bound=1
    )
    run_all_simulations_for_multiple_parent_directory_with_tree_pedigree_subdirectories(
        parent_directory=parent_directory,
        results_path=results_directory_path,
        merge_error_operations_per_simulation=merge_error_operations_per_simulation
    )


def run_interactive_session():
    running_mode_input = ("Specify the running mode:\n"
                          "1) Specify a single tree-pedigree directory\n"
                          "2) Specify a parent directory containing tree-pedigree subdirectories\n"
                          "3) Specify a super-parent directory whose subdirectories contain tree-pedigree directories "
                          "(Run the previous mode for all the subdirectories of a directory)\n")
    running_mode = get_natural_number_input_in_bounds(running_mode_input, 1, 3)
    match running_mode:
        case 1:
            run_interactive_session_single_directory()
        case 2:
            run_interactive_session_multiple_tree_pedigree_subdirectories_all_trees()
        case 3:
            run_interactive_session_multiple_parent_directories_with_tree_pedigree_subdirectories_all_trees()


if __name__ == '__main__':
    run_interactive_session()
