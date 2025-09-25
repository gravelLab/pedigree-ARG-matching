import os
import shutil
import warnings
from pathlib import Path
from typing import Iterable

from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.configuration import tree_extension
from scripts.simulations.tree_error_simulation.simulate_resolve_polytomy_tree_error import oracle_filename
from scripts.utility.basic_utility import (get_directory_path, get_non_existing_path,
                                           get_unique_filename_with_specified_extension, get_filepath,
                                           get_number_input_in_bounds)
from enum import Enum

metadata_header = ("This file specifies the number of probands for every tree in this directory\n"
                   "tree_index,proband_number")
metadata_filename = ".metadata.txt"

minimum_proband_number = 8


class Direction(Enum):
    UP = 1,
    DOWN = 2


def process_tree(tree_filepath: Path, tree_result_directory: Path):
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
    vertices_to_remove = [x for x in original_coalescent_tree if original_coalescent_tree.has_parents(x)]
    # Process every edge
    for vertex in vertices_to_remove:
        # Keeping both bottom and upper parts of the tree
        trees: Iterable[CoalescentTree] = original_coalescent_tree.subdivide_tree(edge_child_vertex=vertex)
        for tree in trees:
            if len(tree.get_sink_vertices()) >= minimum_proband_number:
                result_trees.append(tree)
    # Save the trees by the number of probands
    result_trees.sort(key=lambda cut_tree: len(cut_tree.get_sink_vertices()), reverse=True)
    metadata_filepath = tree_result_directory / metadata_filename
    with open(metadata_filepath, "w") as metadata_file:
        for i, tree in enumerate(result_trees):
            modified_tree_filepath = tree_result_directory / f"{i}"
            tree.save_to_file(filepath=modified_tree_filepath)
            metadata_file.write(f"{i}: {len(tree.get_sink_vertices())}\n")
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


def process_single_tree():
    tree_path = Path(get_filepath("Specify the path to the tree: "))
    tree_path_parent_dir = tree_path.parent
    os.chdir(tree_path_parent_dir)
    result_directory_name = get_non_existing_path("Specify the resulting directory name/path:")
    result_directory_path = tree_path_parent_dir / result_directory_name
    process_tree(tree_path, result_directory_path)


def process_proband_sorted_trees_for_error_simulations():
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


def run_interactive_session():
    menu_prompt = ("Specify the running mode:\n"
                   "1) Create subclades from a single tree\n"
                   "2) (Tree error simulations) Specify the path to a parent directory with trees grouped by the "
                   "proband number\n")
    menu_option = get_number_input_in_bounds(menu_prompt, 1, 2)
    match menu_option:
        case 1:
            process_single_tree()
        case 2:
            process_proband_sorted_trees_for_error_simulations()


if __name__ == '__main__':
    run_interactive_session()
