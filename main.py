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
from divide_and_conquer.graph_matcher import GraphMather, MatcherLogger, logs_default_directory_name
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from divide_and_conquer.subtree_matcher import SubtreeMatcher


def get_file_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        else:
            return file_path


def get_non_existing_directory_name(input_request: str):
    while True:
        directory_name = input(input_request)
        if os.path.exists(directory_name):
            print("The specified directory already exists, try again")
        else:
            return directory_name


def save_dictionary_to_file(dictionary_filename: str, dictionary: dict):
    dictionary_file = open(dictionary_filename, 'w')
    for key, value in dictionary.items():
        dictionary_file.write(f"{key}: {value}\n")
    dictionary_file.close()


def save_statistics_to_file(alignment_result: {int: [SubtreeMatcher]}, coalescent_tree: CoalescentTree,
                            clade_root: int, alignments_number: int, filename: str):
    # Calculating the clade for the given vertex and sorting the vertices by their levels
    clade = coalescent_tree.get_connected_component_for_vertex(clade_root)
    clade = sorted(clade, key=lambda v: coalescent_tree.vertex_to_level_map[v], reverse=True)
    # Calculating the rest of the data
    coalescing_events = len([x for x in clade if x in coalescent_tree.children_map
                             and len(coalescent_tree.children_map[x]) > 1])
    coalescent_vertex_pedigree_candidates_number = {x: len(alignment_result[x]) for x in clade}
    # Printing the results to the file
    statistics_file = open(filename, 'w')
    statistics_file.write(f"The root of the clade: {clade_root}\n")
    statistics_file.write(f"There are {len(clade)} vertices in the clade\n")
    statistics_file.write(f"There are {coalescing_events} coalescing events in the clade\n")
    statistics_file.write(f"The total number of alignments is: {alignments_number}\n")
    statistics_file.write("##############################\n")
    statistics_file.write("The number of pedigree candidates for every vertex:\n")
    for vertex in clade:
        statistics_file.write(f"{vertex}: {coalescent_vertex_pedigree_candidates_number[vertex]}\n")
    statistics_file.close()


def save_alignment_result_to_files(alignment_result: {int: [SubtreeMatcher]}, coalescent_tree: CoalescentTree):
    top_level_vertices = coalescent_tree.levels[-1]
    for top_level_vertex in top_level_vertices:
        # TODO: Draw an image?
        directory_name = f"{top_level_vertex}"
        os.makedirs(directory_name)
        os.chdir(directory_name)
        counter = 0
        subtree_matchers = alignment_result[top_level_vertex]
        for subtree_matcher in subtree_matchers.values():
            alignments = subtree_matcher.get_all_subtree_alignments()
            for alignment in alignments:
                alignment_filename = f"alignment_{counter}"
                save_dictionary_to_file(dictionary_filename=alignment_filename, dictionary=alignment)
                counter += 1
        statistics_file_name = f"clade_{top_level_vertex}.txt"
        save_statistics_to_file(alignment_result=alignment_result, coalescent_tree=coalescent_tree,
                                clade_root=top_level_vertex, alignments_number=counter, filename=statistics_file_name)
        os.chdir("..")


os.chdir("pedigrees")
directory = get_file_path("Specify the testing directory in the pedigrees folder:")
os.chdir(directory)
result_directory_name = get_non_existing_directory_name("Specify the name of the directory containing the result:")
os.makedirs(result_directory_name)

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
    # Skipping the directories
    if os.path.isdir(file):
        continue
    filename = os.path.basename(file)
    if filename == f"{directory}.pedigree":
        continue
    print(f"{filename}")
    logs_directory_name = f"{result_directory_name}/{filename}/{logs_default_directory_name}"
    os.makedirs(logs_directory_name)
    logger = MatcherLogger(logs_directory_name=logs_directory_name)
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filename)
    graph_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                logger=logger)
    start_alignment = time.time()
    result = graph_matcher.find_mapping()
    end_alignment = time.time()
    print(f"Matching time: {end_alignment - start_alignment} seconds")
    total_alignment_time += end_alignment - start_alignment
    count += 1
    os.chdir(result_directory_name)
    os.chdir(filename)
    save_alignment_result_to_files(result, coalescent_tree)
    os.chdir("../..")
print(f"Average inference time {total_alignment_time / count} seconds")
