"""!
@file main.py
@brief The main script that asks the user to specify a testing directory for the script. The testing directory
must be located under the "pedigrees" folder and the user should only specify the name of the subdirectory.
The testing directory must contain one file with the extension ".pedigree" and the rest of the files must be the
coalescent trees with which the mentioned pedigree is to be aligned.
"""

import os
import time

from genealogical_graph import *
from divide_and_conquer.graph_matcher import GraphMather, MatcherLogger
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph

directory = input("Specify the testing directory in the pedigrees folder:")
os.chdir(f"pedigrees/{directory}")

# filename_input = input("Specify the pedigree file:")
filename = f"{directory}.pedigree"

pedigree = Graph.get_pedigree_from_file(filename=filename)
start_preprocessing = time.time()
potential_mrca_graph = PotentialMrcaProcessedGraph(pedigree=pedigree)
end_preprocessing = time.time()

print(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")

# Get the current working directory
current_directory = os.getcwd()
# Get a list of all files in the current directory
total_alignment_time = 0
count = 0
for file in os.listdir(current_directory):
    filename = os.path.basename(file)
    if filename == f"{directory}.pedigree":
        continue
    print(f"{filename}")
    logger = MatcherLogger()
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filename)
    graph_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                logger=logger)
    start_alignment = time.time()
    graph_matcher.find_mapping()
    end_alignment = time.time()
    print(f"Matching time: {end_alignment - start_alignment} seconds")
    total_alignment_time += end_alignment - start_alignment
    count += 1
print(f"Average inference time {total_alignment_time / count} seconds")
