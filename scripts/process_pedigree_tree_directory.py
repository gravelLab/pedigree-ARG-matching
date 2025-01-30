"""!
@file process_pedigree_tree_directory.py
@brief The main script that asks the user to specify a testing directory for the script.
The testing directory must contain one file with the extension ".pedigree" and the rest of the files must be the
coalescent trees with which the mentioned pedigree is to be aligned.
"""
import os

from scripts.run_alignment import process_pedigree_tree_directory
from scripts.utility import get_directory_path, get_non_existing_path


def run_interactive_session():
    pedigree_directory = get_directory_path("Specify the testing directory:")
    os.chdir(pedigree_directory)
    result_dir_name = get_non_existing_path("Specify the path to the directory to save the result:")
    process_pedigree_tree_directory(pedigree_directory, result_dir_name)


if __name__ == "__main__":
    run_interactive_session()
