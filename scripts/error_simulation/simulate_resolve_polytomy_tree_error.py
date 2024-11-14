from pathlib import Path

from graph.coalescent_tree import CoalescentTree
from scripts.utility import *


def run_independent_simulations(tree: CoalescentTree, simulation_number: int, polytomy_vertices: [int],
                                tree_absolute_path: str, result_directory_path: Path):
    simulation_error_vertices = list(polytomy_vertices)
    for simulation_counter in range(simulation_number):
        simulation_tree_name = str(simulation_counter)
        simulation_tree_result_path = result_directory_path / simulation_tree_name
        error_vertex = random.choice(simulation_error_vertices)
        # Removing the vertex from the list to prevent it from being picked again
        simulation_error_vertices.remove(error_vertex)
        tree.unmerge_polytomy(child=error_vertex, recalculate_levels=True)
        tree.save_to_file(filepath=simulation_tree_result_path)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_absolute_path)


def run_all_simulations(tree: CoalescentTree, polytomy_vertices: [int],
                        tree_absolute_path: str, result_directory_path: Path):
    for simulation_counter, error_vertex in enumerate(polytomy_vertices):
        simulation_tree_name = str(simulation_counter)
        simulation_tree_result_path = result_directory_path / simulation_tree_name
        tree.unmerge_polytomy(child=error_vertex, recalculate_levels=True)
        tree.save_to_file(filepath=simulation_tree_result_path)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_absolute_path)


def run_interactive_session():
    result_directory_paths = Path(get_non_existing_directory_path("Specify the path where the"
                                                                  " results are to be stored:"))
    coalescent_tree_path = get_file_path("Specify the path to the coalescent tree:")
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
    os.makedirs(result_directory_paths)
    polytomy_parents = (
        parent
        for parent, children in coalescent_tree.children_map.items()
        if len(children) > 2
    )

    possible_error_vertices = [
        vertex
        for parent in polytomy_parents
        for vertex in coalescent_tree.children_map[parent]
    ]

    error_mode_input = ("Specify the error-simulation mode:\n"
                        "1) Perform given number of 1 independent polytomy resolution\n"
                        "2) Generate all the possible 1 polytomy resolutions\n")

    error_mode = get_natural_number_input_in_bounds(error_mode_input, 1, 2)
    if error_mode == 1:
        possible_trees_number = len(possible_error_vertices)
        custom_simulation_number_input_request = (f"Specify the number of simulations. Notice that there are "
                                                  f"{possible_trees_number} possible trees with 1"
                                                  f" resolved polytomy\n"
                                                  f"Any number higher than that will be treated by simulating "
                                                  f"all the possible errors (the 2 option)")
        simulation_number = get_natural_number_input(custom_simulation_number_input_request)
        if simulation_number > possible_trees_number:
            run_all_simulations(tree=coalescent_tree, polytomy_vertices=possible_error_vertices,
                                tree_absolute_path=coalescent_tree_path, result_directory_path=result_directory_paths)
        else:
            run_independent_simulations(tree=coalescent_tree, simulation_number=simulation_number,
                                        polytomy_vertices=possible_error_vertices,
                                        tree_absolute_path=coalescent_tree_path,
                                        result_directory_path=result_directory_paths)
    else:
        run_all_simulations(tree=coalescent_tree, polytomy_vertices=possible_error_vertices,
                            tree_absolute_path=coalescent_tree_path, result_directory_path=result_directory_paths)


if __name__ == "__main__":
    run_interactive_session()
