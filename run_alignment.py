import os
import time

from genealogical_graph import *
from divide_and_conquer.graph_matcher import GraphMather, MatcherLogger, logs_default_directory_name
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from utility import *

print_enabled = False


def save_dictionary_to_file(dictionary_filename: str, dictionary: dict):
    dictionary_file = open(dictionary_filename, 'w')
    for key, value in dictionary.items():
        dictionary_file.write(f"{key}: {value}\n")
    dictionary_file.close()


def save_statistics_to_file(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                            clade_root: int, alignments_number: int, filename: str):
    # Calculating the clade for the given vertex and sorting the vertices by their levels
    clade = coalescent_tree.get_connected_component_for_vertex(clade_root)
    clade = sorted(clade, key=lambda v: coalescent_tree.vertex_to_level_map[v], reverse=True)
    # Calculating the rest of the data
    coalescing_events = len([x for x in clade if x in coalescent_tree.children_map
                             and len(coalescent_tree.children_map[x]) > 1])
    proband_vertices = [x for x in clade if coalescent_tree.vertex_to_level_map[x] == 0]
    coalescent_vertex_pedigree_candidates_number = {x: len({alignment[x] for alignment in alignment_result[clade_root]})
                                                    for x in clade}
    # Printing the results to the file
    statistics_file = open(filename, 'w')
    statistics_file.write(f"The root of the clade: {clade_root}\n")
    statistics_file.write(f"There are {len(clade)} vertices in the clade\n")
    statistics_file.write(f"There are {len(proband_vertices)} probands in the clade\n")
    statistics_file.write(f"There are {coalescing_events} coalescing events in the clade\n")
    statistics_file.write("##############################\n")
    statistics_file.write("Number of coalescing events grouped by the children number:\n")
    coalescing_number_to_events_number = dict()
    for vertex in clade:
        if vertex not in coalescent_tree.children_map:
            continue
        children_number = len(coalescent_tree.children_map[vertex])
        previous_counter = coalescing_number_to_events_number.get(children_number, 0)
        coalescing_number_to_events_number[children_number] = previous_counter + 1
    coalescing_numbers = sorted(coalescing_number_to_events_number.keys())
    for coalescing_number in coalescing_numbers:
        statistics_file.write(f"{coalescing_number}: {coalescing_number_to_events_number[coalescing_number]}\n")
    statistics_file.write("##############################\n")
    statistics_file.write(f"The total number of alignments is: {alignments_number}\n")
    statistics_file.write("##############################\n")
    statistics_file.write("The number of pedigree candidates for every vertex:\n")
    for vertex in clade:
        statistics_file.write(f"{vertex}: {coalescent_vertex_pedigree_candidates_number[vertex]}\n")
    statistics_file.close()


def save_general_statistics_coalescent_tree(coalescent_tree: CoalescentTree, coalescent_tree_filename: str):
    statistics_file = open(coalescent_tree_filename, "w")
    top_level_vertices = coalescent_tree.levels[-1]
    statistics_file.write(f"There are {len(top_level_vertices)} clades in the coalescent tree\n")
    statistics_file.close()


def save_alignment_result_to_files(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                                   coalescent_tree_filename: str):
    save_general_statistics_coalescent_tree(coalescent_tree, coalescent_tree_filename)
    for clade_root in alignment_result:
        # TODO: Draw an image?
        directory_name = f"{clade_root}"
        os.makedirs(directory_name)
        os.chdir(directory_name)
        counter = 0
        valid_alignments = alignment_result[clade_root]
        for valid_alignment in valid_alignments:
            alignment_filename = f"alignment_{counter}"
            save_dictionary_to_file(dictionary_filename=alignment_filename, dictionary=valid_alignment)
            counter += 1
        statistics_file_name = f"clade_{clade_root}.txt"
        save_statistics_to_file(alignment_result=alignment_result, coalescent_tree=coalescent_tree,
                                clade_root=clade_root, alignments_number=counter, filename=statistics_file_name)
        os.chdir("..")


def run_alignment_and_save_results(directory: str, result_directory_name: str):
    os.chdir(directory)
    os.makedirs(result_directory_name)
    # filename_input = input("Specify the pedigree file:")
    pedigree_files = [file for file in os.listdir() if file.endswith('.pedigree')]
    if len(pedigree_files) == 0:
        raise Exception("There are no pedigree files in the directory")
    if len(pedigree_files) == 2:
        raise Exception("There are multiple pedigree files in the directory")

    pedigree_filename = pedigree_files[0]

    pedigree = Graph.get_pedigree_from_file(filename=pedigree_filename)
    start_preprocessing = time.time()
    potential_mrca_graph = PotentialMrcaProcessedGraph(pedigree=pedigree)
    end_preprocessing = time.time()
    if print_enabled:
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
        if filename == pedigree_filename:
            continue
        if print_enabled:
            print(f"{filename}")
        logs_directory_name = f"{result_directory_name}/{filename}/{logs_default_directory_name}"
        os.makedirs(logs_directory_name)
        logger = MatcherLogger(logs_directory_name=logs_directory_name)
        coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filename)
        coalescent_tree.remove_unary_nodes()
        graph_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                    logger=logger)
        start_alignment = time.time()
        result = graph_matcher.find_mapping()
        end_alignment = time.time()
        if print_enabled:
            print(f"Matching time: {end_alignment - start_alignment} seconds")
        total_alignment_time += end_alignment - start_alignment
        count += 1
        os.chdir(result_directory_name)
        os.chdir(filename)
        save_alignment_result_to_files(result, coalescent_tree, filename)
        os.chdir("../..")
    if print_enabled:
        print(f"Average inference time {total_alignment_time / count} seconds")
    os.chdir("..")
