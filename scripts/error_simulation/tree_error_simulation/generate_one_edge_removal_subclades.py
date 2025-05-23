import os
import shutil
import warnings
from pathlib import Path

from alignment.configuration import tree_extension
from graph.coalescent_tree import CoalescentTree, Direction
from scripts.error_simulation.tree_error_simulation.simulate_resolve_polytomy_tree_error import oracle_filename
from scripts.utility.basic_utility import (get_directory_path, get_non_existing_path,
                                           get_unique_filename_with_specified_extension)

metadata_header = ("This file specifies the number of probands for every tree in this directory\n"
                   "tree_index,proband_number")
metadata_filename = ".metadata.txt"

minimum_proband_number = 8


def process_tree(tree_filepath: Path, tree_result_directory):
    try:
        original_coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_filepath)
    except Exception as ex:
        print(f"Could not process {tree_filepath} due to an exception:")
        print(ex)
        return
    tree_parent_directory = tree_filepath.parent
    os.makedirs(tree_result_directory, exist_ok=True)
    # Saving the original tree as well
    result_trees = [original_coalescent_tree]
    # Finding the edges to be cut. Notice that in a tree an edge can be identified by the child id
    vertices_to_remove = [x for x in original_coalescent_tree.parents_map if
                          original_coalescent_tree.parents_map.get(x, [])]
    # Process every edge
    for vertex in vertices_to_remove:
        # Keeping both bottom and upper parts of the tree
        for direction in (Direction.UP, Direction.DOWN):
            cloned_coalescent_tree = original_coalescent_tree.clone()
            cloned_coalescent_tree.cut_edge(edge_child_vertex=vertex, direction=direction)
            if len(cloned_coalescent_tree.probands) >= minimum_proband_number:
                result_trees.append(cloned_coalescent_tree)
    # Save the trees by the number of probands
    result_trees.sort(key=lambda cut_tree: len(cut_tree.get_probands()), reverse=True)
    metadata_filepath = tree_result_directory / metadata_filename
    with open(metadata_filepath, "w") as metadata_file:
        for i, tree in enumerate(result_trees):
            modified_tree_filepath = tree_result_directory / f"{i}"
            tree.save_to_file(filepath=modified_tree_filepath)
            metadata_file.write(f"{i}: {len(tree.probands)}\n")
    oracle_origin_filepath = tree_parent_directory / oracle_filename
    if os.path.exists(oracle_origin_filepath):
        oracle_result_filepath = tree_result_directory / tree_result_directory
        shutil.copy(oracle_origin_filepath, oracle_result_filepath)


def process_tree_directory(tree_directory, tree_directory_result_filepath):
    for file in tree_directory.iterdir():
        if not file.is_dir():
            continue
        try:
            tree_filename = get_unique_filename_with_specified_extension(directory_path=file, extension=tree_extension)
        except ValueError as _:
            warnings.warn(f"Could not find the tree file under {tree_directory}")
            continue
        tree_filepath = file / tree_filename
        trees_result_directory = tree_directory_result_filepath / file.name
        process_tree(tree_filepath=tree_filepath, tree_result_directory=trees_result_directory)


def process_proband_directory(proband_directory, proband_result_filepath):
    for tree_directory in proband_directory.iterdir():
        if not tree_directory.is_dir():
            # Skipping all the files
            continue
        tree_dir_path = proband_directory / tree_directory.name
        result_tree_dir_path = proband_result_filepath / tree_directory.name
        process_tree_directory(tree_dir_path, result_tree_dir_path)


def run_interactive_session():
    parent_dir_filepath = get_directory_path("Specify the path to the parent directory with error simulated trees:")
    parent_dir_filepath = Path(parent_dir_filepath)
    result_path = get_non_existing_path("Specify the resulting path:")
    result_path = Path(result_path)
    os.makedirs(result_path, exist_ok=True)
    for proband_directory in os.listdir(parent_dir_filepath):
        proband_filepath = parent_dir_filepath / proband_directory
        proband_result_filepath = result_path / proband_directory
        if not os.path.isdir(proband_filepath):
            continue
        try:
            # noinspection PyUnusedLocal
            proband_number = int(proband_directory)
        except ValueError:
            continue
        process_proband_directory(proband_directory=proband_filepath,
                                  proband_result_filepath=proband_result_filepath)


if __name__ == '__main__':
    run_interactive_session()
