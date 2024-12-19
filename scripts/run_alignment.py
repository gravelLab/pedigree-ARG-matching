import sys
from itertools import combinations

from alignment.graph_matcher import *
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.genealogical_graph import *
from scripts.utility import *

print_enabled = False
calculate_similarity = False
pedigree_file_extension = ".pedigree"
extension_skip_list = [".svg", pedigree_file_extension]


def get_assignment_likelihood(coalescent_tree: CoalescentTree, alignment: dict, pedigree: SimpleGraph):
    def get_coalescent_vertex_mapped_level_difference(ancestor: int, descendant: int) -> int:
        return pedigree.get_minimal_path_length(ancestor=alignment[ancestor], descendant=alignment[descendant])

    total_length = 0
    for parent, children in coalescent_tree.children_map.items():
        if parent not in alignment:
            # The vertex does not belong to this clade
            continue
        for child in children:
            level_distance_approximation = get_coalescent_vertex_mapped_level_difference(ancestor=parent,
                                                                                         descendant=child)
            total_length += level_distance_approximation
    return total_length


def save_alignment_result_to_file(coalescent_tree: CoalescentTree, pedigree: GenealogicalGraph,
                                  dictionary_filepath: str, alignment: dict):
    alignment_likelihood = get_assignment_likelihood(coalescent_tree, alignment, pedigree)
    dictionary_file = open(dictionary_filepath, 'w')
    dictionary_file.write(f"Approximated alignment length: {alignment_likelihood}\n")
    for key, value in alignment.items():
        dictionary_file.write(f"{key}: {value}\n")
    dictionary_file.close()
    return alignment_likelihood


def get_alignments_proband_distance_probands_ignored(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    # Calculate the number of shared values. Assuming here that the keys in both the dictionaries are the same
    distance = sum(1 for key in first if key not in probands and first[key] != second[key])
    assert distance == len(first) - len(probands) - sum(1 for key in first if key not in probands and
                                                        first[key] == second[key])
    return distance


def get_alignments_key_similarity_probands_ignored(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    # Calculate the number of shared values. Assuming here that the keys in both the dictionaries are the same
    return sum(1 for key in first if key not in probands and first[key] == second[key])


def get_alignments_pair_similarity(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    probands = set(probands).intersection(first)
    # Calculate the number of shared values
    shared_values = get_alignments_key_similarity_probands_ignored(first, second, probands)

    # Calculate the total size (number of elements)
    total_size = len(first) - len(probands)

    # Calculate the similarity using the formula
    similarity = shared_values / total_size

    return similarity


def get_alignments_similarity(alignments: [dict], probands):
    total_similarity = 0
    for dict1, dict2 in combinations(alignments, 2):
        total_similarity += get_alignments_pair_similarity(first=dict1, second=dict2, probands=probands)
    return total_similarity / (len(alignments) * (len(alignments) - 1) / 2)


def get_alignments_pair_ind_similarity(first: dict, second: dict, probands: {int}):
    # Ensure both dictionaries have the same number of elements
    if len(first) != len(second):
        raise ValueError("Both dictionaries must have the same number of elements.")
    probands = set(probands).intersection(first)
    # Calculate the number of shared values
    shared_values = sum(1 for key in first if key not in probands and first[key] // 2 == second[key] // 2)

    # Calculate the total size (number of elements)
    total_size = len(first) - len(probands)

    # Calculate the similarity using the formula
    similarity = shared_values / total_size

    return similarity


def get_alignments_ind_similarity(alignments: [dict], probands):
    total_similarity = 0
    for dict1, dict2 in combinations(alignments, 2):
        total_similarity += get_alignments_pair_ind_similarity(first=dict1, second=dict2, probands=probands)
    return total_similarity / (len(alignments) * (len(alignments) - 1) / 2)


def get_distance_histogram_to_identity(alignments: [dict]):
    # Assuming that all the alignments have the same keys
    coalescent_vertices = alignments[0].keys()
    vertices_length = len(coalescent_vertices)
    identity_solution = {x: x for x in coalescent_vertices}
    distance_histogram = defaultdict(int)
    for alignment in alignments:
        alignment_distance = (vertices_length -
                              sum(1 for key in identity_solution if identity_solution[key] == alignment[key]))
        # alignment_distance = get_alignments_pair_distance_probands_ignored(alignment, identity_solution, probands)
        distance_histogram[alignment_distance] += 1
    # Verifying that the identity alignment is always present
    # if current_matching_mode == MatchingMode.ALL_ALIGNMENTS:
    #     assert 0 in distance_histogram
    return distance_histogram


def save_statistics_to_file(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                            clade_root: int, alignments_number: int, result_filepath: str,
                            min_alignment_length: int, min_length_alignments: [int]):
    # Calculating the clade for the given vertex and sorting the vertices by their levels
    clade = coalescent_tree.get_vertex_descendants(clade_root)
    # clade = coalescent_tree.get_connected_component_for_vertex(clade_root)
    clade = sorted(clade, key=lambda v: coalescent_tree.vertex_to_level_map[v], reverse=True)
    # Calculating the rest of the data
    coalescing_events = len([x for x in clade if x in coalescent_tree.children_map
                             and len(coalescent_tree.children_map[x]) > 1])
    proband_vertices = [x for x in clade if coalescent_tree.vertex_to_level_map[x] == 0]
    coalescent_vertex_pedigree_candidates_number = {x: len({alignment[x] for alignment in alignment_result[clade_root]})
                                                    for x in clade}
    # Printing the results to the file
    statistics_file = open(result_filepath, 'w')
    statistics_file.write(f"The root of the clade: {clade_root}\n")
    statistics_file.write(f"There are {len(clade)} vertices in the clade\n")
    statistics_file.write(f"There are {len(proband_vertices)} probands in the clade\n")
    statistics_file.write(f"There are {coalescing_events} coalescing events in the clade\n")
    statistics_file.write("##############################\n")
    statistics_file.write(f"The minimum alignment length: {min_alignment_length}\n")
    statistics_file.write(f"The corresponding alignments: {min_length_alignments}\n")
    statistics_file.write("##############################\n")
    if calculate_similarity:
        alignments_similarity = get_alignments_similarity(alignment_result[clade_root], coalescent_tree.probands)
        alignments_ind_similarity = get_alignments_ind_similarity(alignment_result[clade_root],
                                                                  coalescent_tree.probands)
        statistics_file.write(f"The alignments similarity: {alignments_similarity}\n")
        statistics_file.write(f"The alignments individual similarity: {alignments_ind_similarity}\n")
    else:
        statistics_file.write("The similarities among the trees were not calculated\n")
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
        if vertex in coalescent_tree.probands:
            statistics_file.write(f"{vertex} (proband): {coalescent_vertex_pedigree_candidates_number[vertex]}\n")
        else:
            statistics_file.write(f"{vertex}: {coalescent_vertex_pedigree_candidates_number[vertex]}\n")
    statistics_file.write("##############################\n")
    statistics_file.write("Appearances of each candidate:\n")
    vertex_candidate_count_dict = defaultdict(lambda: defaultdict(int))
    for alignment in alignment_result[clade_root]:
        for key, value in alignment.items():
            vertex_candidate_count_dict[key][value] += 1
    for vertex in clade:
        statistics_file.write(f"{vertex}:\n")
        for value, count in sorted(vertex_candidate_count_dict[vertex].items(), key=lambda x: x[1], reverse=True):
            percentage = count / alignments_number
            statistics_file.write(f"        {value} ({percentage});\n")
    statistics_file.close()


def save_alignment_result_to_files(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                                   pedigree: GenealogicalGraph, directory_path: str = ""):
    result_parent_directory = Path(directory_path)
    for clade_root in alignment_result:
        # TODO: Draw an image?
        result_directory_name = f"{clade_root}"
        clade_result_directory_path = result_parent_directory / result_directory_name
        os.makedirs(clade_result_directory_path)
        # os.chdir(result_directory_name)
        counter = 0
        valid_alignments = alignment_result[clade_root]
        min_length = sys.maxsize
        min_length_alignments = []

        for valid_alignment in valid_alignments:
            alignment_filename = f"alignment_{counter}"
            alignment_path = clade_result_directory_path / alignment_filename
            alignment_length = save_alignment_result_to_file(dictionary_filepath=alignment_path,
                                                             alignment=valid_alignment,
                                                             pedigree=pedigree, coalescent_tree=coalescent_tree)
            if min_length == alignment_length:
                min_length_alignments.append(counter)
            elif alignment_length < min_length:
                min_length = alignment_length
                min_length_alignments = [counter]
            counter += 1
        statistics_filename = f"_clade_{clade_root}.txt"
        statistics_filepath = clade_result_directory_path / statistics_filename
        save_statistics_to_file(alignment_result=alignment_result, coalescent_tree=coalescent_tree,
                                clade_root=clade_root, alignments_number=counter, result_filepath=statistics_filepath,
                                min_alignment_length=min_length, min_length_alignments=min_length_alignments)
        if len(valid_alignments) > 0:
            distance_histogram_filename = f"distance_histogram_{clade_root}.txt"
            distance_histogram_filepath = clade_result_directory_path / distance_histogram_filename
            distance_histogram_dict = get_distance_histogram_to_identity(valid_alignments)
            save_dictionary_to_file(dictionary_filepath=distance_histogram_filepath, dictionary=distance_histogram_dict)
            distance_histogram_image_filename = f"distance_histogram_{clade_root}"
            distance_histogram_image_filepath = clade_result_directory_path / distance_histogram_image_filename
            build_histogram(histogram_filepath=distance_histogram_image_filepath, dictionary=distance_histogram_dict)
        # os.chdir("..")


def process_pedigree_tree_directory(directory: str, result_directory_name: str):
    parent_directory_path = Path(directory)
    result_directory_path = parent_directory_path / result_directory_name
    os.makedirs(result_directory_path)
    directory_files = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]
    if len(directory_files) != 2:
        raise ValueError(f"The pedigree-tree directory must contain only two files, that is the pedigree and "
                         f"the coalescent tree. The pedigree file must have the {pedigree_file_extension} extension")
    pedigree_filename = get_unique_filename_with_specified_extension(directory_path=directory,
                                                                     extension=pedigree_file_extension)
    tree_filename = [filename for filename in directory_files if filename != pedigree_filename][0]
    tree_filepath = parent_directory_path / tree_filename
    start_preprocessing = time.time()
    print_log(f"Processing the coalescent tree: {tree_filepath}")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_filepath)
    coalescent_tree.remove_unary_nodes()
    pedigree_filepath = parent_directory_path / pedigree_filename
    print_log(f"Processing the pedigree: {pedigree_filepath}")
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_filepath,
                                                                         preprocess_graph=False)
    if default_initial_matching_mode == InitialMatchingMode.PLOID:
        pedigree_probands = coalescent_tree.get_probands()
    else:
        proband_individuals = {proband // 2 for proband in coalescent_tree.get_probands()}
        probands_all_ploids = [
            ploid
            for proband_individual in proband_individuals
            for ploid in [2 * proband_individual, 2 * proband_individual + 1]
        ]
        pedigree_probands = probands_all_ploids
    pedigree.reduce_to_ascending_genealogy(probands=pedigree_probands,
                                           recalculate_levels=True)
    pedigree.initialize_potential_mrca_map()
    end_preprocessing = time.time()
    print_log(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")
    alignment_result_path = result_directory_path / tree_filename
    logs_directory_path = alignment_result_path / logs_default_directory_name
    os.makedirs(logs_directory_path)
    logger = MatcherLogger(logs_directory_path=logs_directory_path)
    graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                 logger=logger)
    start_alignment = time.time()
    print_log("Running the alignment")
    result = graph_matcher.find_mapping()
    end_alignment = time.time()
    print_log(f"Matching time: {end_alignment - start_alignment} seconds")
    logger.log(f"Matching time: {end_alignment - start_alignment} seconds")
    save_alignment_result_to_files(alignment_result=result, coalescent_tree=coalescent_tree,
                                   pedigree=pedigree, directory_path=str(alignment_result_path))


def run_alignment_with_multiple_clades_and_save_results(directory: str, result_directory_name: str,
                                                        pedigree: PotentialMrcaProcessedGraph = None):
    parent_directory_path = Path(directory)
    result_directory_path = parent_directory_path / result_directory_name
    os.makedirs(result_directory_path)
    if pedigree is None:
        pedigree_filename = get_unique_filename_with_specified_extension(directory_path=directory,
                                                                         extension=pedigree_file_extension)
        pedigree_filepath = parent_directory_path / pedigree_filename
        start_preprocessing = time.time()
        pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_filepath)
        end_preprocessing = time.time()
        print_log(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")
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
        extension = os.path.splitext(filename)[1]
        if extension in extension_skip_list or filename == "SKIP":
            continue
        print_log(f"{filename}")
        alignment_result_path = result_directory_path / filename
        logs_directory_path = alignment_result_path / logs_default_directory_name
        os.makedirs(logs_directory_path)
        logger = MatcherLogger(logs_directory_path=logs_directory_path)
        absolute_path = os.path.abspath(file)
        coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=absolute_path)
        coalescent_tree.remove_unary_nodes()
        graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                     logger=logger)
        start_alignment = time.time()
        result = graph_matcher.find_mapping()
        end_alignment = time.time()
        print_log(f"Matching time: {end_alignment - start_alignment} seconds")
        logger.log(f"Matching time: {end_alignment - start_alignment} seconds")
        total_alignment_time += end_alignment - start_alignment
        count += 1
        save_alignment_result_to_files(alignment_result=result, coalescent_tree=coalescent_tree,
                                       pedigree=pedigree, directory_path=str(alignment_result_path))
    print_log(f"Average inference time {total_alignment_time / count} seconds")
