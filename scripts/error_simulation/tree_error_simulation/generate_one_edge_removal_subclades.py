import os
from pathlib import Path

from graph.coalescent_tree import CoalescentTree
from scripts.utility import get_directory_path, get_non_existing_path

metadata_header = ("This file specifies the number of probands for every tree in this directory\n"
                   "tree_index,proband_number")
metadata_filename = ".metadata.txt"


def process_tree(filepath, tree_result_directory):
    try:
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    except Exception as ex:
        print(f"Could not process {filepath} due to an exception")
        print(ex)
        return
    os.makedirs(tree_result_directory, exist_ok=True)
    # Saving the original tree as well
    result_trees = [coalescent_tree]
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    # Finding the edges to be cut
    vertices_to_remove = [x for x in coalescent_tree.parents_map if coalescent_tree.parents_map.get(x, [])]
    # Process every edge
    for vertex in vertices_to_remove:
        coalescent_tree.cut_edge(edge_child_vertex=vertex)
        if coalescent_tree.probands:
            result_trees.append(coalescent_tree)
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    # Save the trees by the number of probands
    result_trees.sort(key=lambda cut_tree: len(cut_tree.get_probands()), reverse=True)
    metadata_filepath = tree_result_directory / metadata_filename
    with open(metadata_filepath, "w") as metadata_file:
        for i, tree in enumerate(result_trees):
            tree_filepath = tree_result_directory / f"{i}"
            tree.save_to_file(filepath=tree_filepath)
            metadata_file.write(f"{i}: {len(tree.probands)}\n")


def process_tree_directory(tree_directory, tree_directory_result_filepath):
    for file in tree_directory.iterdir():
        if file.is_dir():
            continue
        filepath = file.absolute()
        trees_result_directory = tree_directory_result_filepath / file.name
        process_tree(filepath, trees_result_directory)


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
            proband_number = int(proband_directory)
        except ValueError:
            continue
        process_proband_directory(proband_directory=proband_filepath,
                                  proband_result_filepath=proband_result_filepath)


if __name__ == '__main__':
    run_interactive_session()
