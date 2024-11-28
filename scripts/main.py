"""!
@file main.py
@brief The main script that asks the user to specify a testing directory for the script.
The testing directory must contain one file with the extension ".pedigree" and the rest of the files must be the
coalescent trees with which the mentioned pedigree is to be aligned.
"""
from scripts.run_alignment import run_alignment_and_save_results
from scripts.utility import get_directory_path, get_non_existing_path
import os

pedigree_directory = get_directory_path("Specify the testing directory:")
os.chdir(pedigree_directory)
result_dir_name = get_non_existing_path("Specify the path to the directory to save the result:")
run_alignment_and_save_results(pedigree_directory, result_dir_name)
