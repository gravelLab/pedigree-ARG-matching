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


def run_interactive_single_tree_mode(result_directory_paths: Path):
    coalescent_tree_path = get_filepath("Specify the path to the coalescent tree:")
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
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
    if not possible_error_vertices:
        print("This coalescent tree has not polytomies to resolve")
        return
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


def run_interactive_multiple_parent_tree_directories(result_directory_paths: Path):
    parent_trees_directory = Path(get_directory_path("Specify the path to the parent trees directory:"))

    pass


def run_interactive_session():
    result_directory_paths = Path(get_non_empty_string("Specify the path where the results are to be stored:"))
    os.makedirs(result_directory_paths, exist_ok=True)
    mode_prompt = """
                  Specify the running mode:
                  1) Run the error simulation for a single tree.
                  2) Specify the path to a parent directory containing the tree-pedigree subdirectories. The error 
                  simulation will be performed for all the trees.
                  """
    mode = get_natural_number_input_in_bounds(input_request=mode_prompt,
                                              lower_bound=1,
                                              upper_bound=2)
    match mode:
        case 1:
            run_interactive_single_tree_mode(result_directory_paths=result_directory_paths)
        case 2:
            run_interactive_multiple_parent_tree_directories(result_directory_paths=result_directory_paths)


if __name__ == "__main__":
    run_interactive_session()
