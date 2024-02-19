"""!
@file main.py
@brief The main script that asks the user to specify a testing directory for the script. The testing directory
must be located under the "pedigrees" folder and the user should only specify the name of the subdirectory.
The testing directory must contain one file with the extension ".pedigree" and the rest of the files must be the
coalescent trees with which the mentioned pedigree is to be aligned.
"""
from run_alignment import run_alignment_and_save_results
from utility import *

os.chdir("pedigrees")
pedigree_directory = get_directory_path("Specify the testing directory in the pedigrees folder:")
result_dir_name = get_non_existing_directory_name("Specify the name of the directory containing the result:")
run_alignment_and_save_results(pedigree_directory, result_dir_name)
