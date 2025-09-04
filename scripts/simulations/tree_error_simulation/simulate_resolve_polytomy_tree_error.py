import copy

from lineagekit.core.CoalescentTree import CoalescentTree

from scripts.utility.basic_utility import *
from scripts.utility.alignment_utility import save_dictionary_to_file

oracle_filename = "oracle.txt"
tree_filename = "tree.tree"


def unmerge_polytomy_and_save_results(tree: CoalescentTree, error_vertex: int, simulation_dir_path: str | Path,
                                      oracle_vertex_assignments: dict):
    os.makedirs(simulation_dir_path, exist_ok=True)
    simulation_tree_result_path = simulation_dir_path / tree_filename
    oracle_filepath = simulation_dir_path / oracle_filename
    error_vertex_parent = tree.get_vertex_parent(error_vertex)
    new_vertex_id = tree.unmerge_edge(child=error_vertex)
    simulation_oracle_vertex_assignments = copy.deepcopy(oracle_vertex_assignments)
    simulation_oracle_vertex_assignments[new_vertex_id] = oracle_vertex_assignments[error_vertex_parent]
    tree.save_to_file(filepath=simulation_tree_result_path)
    # For various simulations with msprime data, we need to know the corresponding pedigree vertex
    # assignments for the tree vertices. Considering that resolving a polytomy requires adding a new vertex to the
    # tree, we need to save this information somewhere
    save_dictionary_to_file(dictionary_filepath=oracle_filepath, dictionary=simulation_oracle_vertex_assignments)


def run_independent_simulations(tree: CoalescentTree, simulation_number: int, polytomy_vertices: [int],
                                tree_absolute_path: str, result_directory_path: Path):
    tree_vertices = tree.nodes
    oracle_vertex_assignments = {x: x for x in tree_vertices}
    simulation_error_vertices = list(polytomy_vertices)
    for simulation_counter in range(simulation_number):
        simulation_tree_name = str(simulation_counter)
        simulation_tree_directory = result_directory_path / simulation_tree_name
        os.makedirs(simulation_tree_directory, exist_ok=True)
        error_vertex = random.choice(simulation_error_vertices)
        # Removing the vertex from the list to prevent it from being picked again
        simulation_error_vertices.remove(error_vertex)
        unmerge_polytomy_and_save_results(tree=tree, error_vertex=error_vertex,
                                          simulation_dir_path=simulation_tree_directory,
                                          oracle_vertex_assignments=oracle_vertex_assignments)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_absolute_path)


def run_all_simulations(tree: CoalescentTree, polytomy_vertices: [int],
                        tree_absolute_path: str, result_directory_path: Path):
    tree_vertices = tree.nodes
    oracle_vertex_assignments = {x: x for x in tree_vertices}
    for simulation_counter, error_vertex in enumerate(polytomy_vertices):
        simulation_tree_name = str(simulation_counter)
        simulation_tree_directory = result_directory_path / simulation_tree_name
        unmerge_polytomy_and_save_results(tree=tree, error_vertex=error_vertex,
                                          simulation_dir_path=simulation_tree_directory,
                                          oracle_vertex_assignments=oracle_vertex_assignments)
        # TODO: Reverse the changes instead of parsing the tree again
        tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_absolute_path)


def run_interactive_single_tree_mode(results_directory_paths: Path):
    coalescent_tree_path = get_filepath("Specify the path to the coalescent tree:")
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
    polytomy_parents = (
        node
        for node, parents_dict in coalescent_tree.pred.items()
        if len(parents_dict) > 2
    )

    possible_error_vertices = [
        vertex
        for parent in polytomy_parents
        for vertex in coalescent_tree.get_children(parent)
    ]
    if not possible_error_vertices:
        print("This coalescent tree has not polytomies to resolve")
        return
    error_mode_input = ("Specify the error-simulation mode:\n"
                        "1) Perform given number of 1 independent polytomy resolution\n"
                        "2) Generate all the possible 1 polytomy resolutions\n")

    error_mode = get_number_input_in_bounds(error_mode_input, 1, 2)
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
                                tree_absolute_path=coalescent_tree_path, result_directory_path=results_directory_paths)
        else:
            run_independent_simulations(tree=coalescent_tree, simulation_number=simulation_number,
                                        polytomy_vertices=possible_error_vertices,
                                        tree_absolute_path=coalescent_tree_path,
                                        result_directory_path=results_directory_paths)
    else:
        run_all_simulations(tree=coalescent_tree, polytomy_vertices=possible_error_vertices,
                            tree_absolute_path=coalescent_tree_path, result_directory_path=results_directory_paths)


def run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(parent_directory: Path,
                                                                               results_path: Path):
    for tree_pedigree_directory in os.listdir(parent_directory):
        tree_pedigree_directory_path = parent_directory / tree_pedigree_directory
        if not os.path.isdir(tree_pedigree_directory_path):
            continue
        paths = get_paths_from_tree_pedigree_directory(tree_pedigree_directory_path)
        if not paths:
            continue
        _, tree_path = paths
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path)
        polytomy_parents = (
            node
            for node, children_dict in coalescent_tree.succ.items()
            if len(children_dict) > 2
        )
        possible_error_vertices = [
            vertex
            for parent in polytomy_parents
            for vertex in coalescent_tree.get_children(parent)
        ]
        results_subpath = results_path / tree_pedigree_directory
        os.makedirs(results_subpath, exist_ok=True)
        run_all_simulations(tree=coalescent_tree, tree_absolute_path=tree_path,
                            polytomy_vertices=possible_error_vertices,
                            result_directory_path=results_subpath)


def run_interactive_multiple_parent_tree_directories_all_trees(results_directory_paths: Path):
    parent_trees_directory = Path(get_directory_path("Specify a parent directory containing"
                                                     " tree-pedigree directories:"))
    parent_trees_directory = Path(parent_trees_directory)
    run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(parent_directory=parent_trees_directory,
                                                                               results_path=results_directory_paths)


def run_all_simulations_for_multiple_parent_directory_with_tree_pedigree_subdirectories(parent_directory: Path,
                                                                                        results_directory_paths: Path):
    for tree_pedigrees_parent_directory in os.listdir(parent_directory):
        results_subpath = results_directory_paths / tree_pedigrees_parent_directory
        parent_subdirectory = parent_directory / tree_pedigrees_parent_directory
        if not os.path.isdir(parent_subdirectory):
            continue
        run_all_simulations_for_parent_directory_with_tree_pedigree_subdirectories(
            parent_directory=parent_subdirectory,
            results_path=results_subpath
        )


def run_interactive_session_multiple_parent_directories_with_tree_pedigree_subdirectories_all_trees(
        results_directory_paths: Path):
    parent_directory = get_directory_path("Specify the path to the parent directory with subdirectories containing"
                                          " tree-pedigree directories:")
    parent_directory = Path(parent_directory)
    run_all_simulations_for_multiple_parent_directory_with_tree_pedigree_subdirectories(
        parent_directory=parent_directory,
        results_directory_paths=results_directory_paths
    )


def run_interactive_session():
    mode_prompt = ("Specify the running mode:\n"
                   "1) Run the error simulation for a single tree.\n"
                   "2) Specify the path to a parent directory containing the tree-pedigree subdirectories.\nThe error "
                   "simulation will be performed for all the trees.\n"
                   "3) Specify a directory whose subdirectories contain tree-pedigree directories "
                   "(Run the previous mode for all the subdirectories of a directory).\n")

    mode = get_number_input_in_bounds(input_request=mode_prompt,
                                      lower_bound=1,
                                      upper_bound=3)
    result_directory_paths = Path(get_non_empty_string("Specify the path where the results are to be stored:"))
    os.makedirs(result_directory_paths, exist_ok=True)
    match mode:
        case 1:
            run_interactive_single_tree_mode(results_directory_paths=result_directory_paths)
        case 2:
            run_interactive_multiple_parent_tree_directories_all_trees(results_directory_paths=result_directory_paths)
        case 3:
            run_interactive_session_multiple_parent_directories_with_tree_pedigree_subdirectories_all_trees(
                results_directory_paths=result_directory_paths
            )


if __name__ == "__main__":
    run_interactive_session()
