import sys

from alignment.driver_file import PloidType
from alignment.graph_matcher import *
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.genealogical_graph import *
from scripts.alignment_similarity import get_alignments_similarity, get_alignments_ind_similarity, \
    get_distance_histogram_to_identity, save_distance_histogram
from scripts.utility import *

print_enabled = False
calculate_similarity = False
pedigree_file_extension = ".pedigree"
extension_skip_list = [".svg", pedigree_file_extension]


def get_alignment_likelihood(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                             alignment: dict):
    def get_coalescent_vertex_mapped_level_difference(ancestor: int, descendant: int) -> int:
        return pedigree.get_minimal_path_length(ancestor=alignment[ancestor],
                                                descendant=alignment[descendant]
                                                )

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


def convert_ploid_id_to_individual(ploid_id: int):
    individual_id = ploid_id // 2
    ploid_type = PloidType.Paternal.value if ploid_id % 2 == 0 else PloidType.Maternal.value
    return f"{individual_id}{ploid_type}"


def save_alignment_result_to_file(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                                  dictionary_filepath: str, alignment: dict):
    with open(dictionary_filepath, "w") as dictionary_file:
        alignment_likelihood = None
        if calculate_likelihood:
            alignment_likelihood = get_alignment_likelihood(pedigree=pedigree,
                                                            coalescent_tree=coalescent_tree,
                                                            alignment=alignment)

            dictionary_file.write(f"Approximated alignment length: {alignment_likelihood}\n")
        for key, value in alignment.items():
            converted_value = convert_ploid_id_to_individual(value)
            dictionary_file.write(f"{key}: {converted_value}\n")
        dictionary_file.close()
        return alignment_likelihood


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
    with open(result_filepath, 'w') as statistics_file:
        statistics_file.write(f"The root of the clade: {clade_root}\n")
        statistics_file.write(f"There are {len(clade)} vertices in the clade\n")
        statistics_file.write(f"There are {len(proband_vertices)} probands in the clade\n")
        statistics_file.write(f"There are {coalescing_events} coalescing events in the clade\n")
        statistics_file.write("##############################\n")
        if calculate_distances_histogram:
            statistics_file.write(f"The minimum alignment length: {min_alignment_length}\n")
            statistics_file.write(f"The corresponding alignments: {min_length_alignments}\n")
        statistics_file.write("##############################\n")
        if calculate_similarity:
            alignments_similarity = get_alignments_similarity(alignment_result[clade_root],
                                                              coalescent_tree.probands)
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
            for value, count in sorted(vertex_candidate_count_dict[vertex].items(),
                                       key=lambda x: x[1], reverse=True):
                percentage = count / alignments_number
                statistics_file.write(f"        {value} ({percentage});\n")


def save_alignment_result_to_files(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                                   pedigree: PotentialMrcaProcessedGraph, directory_path: str | Path = ""):
    result_parent_directory = Path(directory_path)
    for clade_root in alignment_result:
        result_directory_name = f"{clade_root}"
        clade_result_directory_path = result_parent_directory / result_directory_name
        os.makedirs(clade_result_directory_path)
        counter = 0
        valid_alignments = alignment_result[clade_root]
        min_length = sys.maxsize
        min_length_alignments = []
        # TODO: Cache the distance between two vertices per pedigree
        for valid_alignment in valid_alignments:
            alignment_filename = f"alignment_{counter}"
            alignment_path = clade_result_directory_path / alignment_filename
            alignment_length = save_alignment_result_to_file(dictionary_filepath=alignment_path,
                                                             alignment=valid_alignment,
                                                             pedigree=pedigree, coalescent_tree=coalescent_tree)
            if calculate_distances_histogram:
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
        if calculate_distances_histogram and len(valid_alignments) > 0:
            distance_histogram_filename = f"distance_histogram_{clade_root}.txt"
            distance_histogram_filepath = clade_result_directory_path / distance_histogram_filename
            distance_histogram_dict = get_distance_histogram_to_identity(valid_alignments)
            save_dictionary_to_file(dictionary_filepath=distance_histogram_filepath,
                                    dictionary=distance_histogram_dict)
            distance_histogram_image_filename = f"distance_histogram_{clade_root}"
            distance_histogram_image_filepath = clade_result_directory_path / distance_histogram_image_filename
            save_distance_histogram(histogram_filepath=distance_histogram_image_filepath,
                                    dictionary=distance_histogram_dict)


def process_pedigree_tree_directory(directory: str, result_directory_name: str):
    parent_directory_path = Path(directory)
    result_directory_path = parent_directory_path / result_directory_name
    os.makedirs(result_directory_path)
    directory_files = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]
    pedigree_filename = get_unique_filename_with_specified_extension(directory_path=directory,
                                                                     extension=pedigree_file_extension)
    tree_filenames = [filename for filename in directory_files if filename != pedigree_filename]
    for tree_filename in tree_filenames:
        tree_filepath = parent_directory_path / tree_filename
        start_preprocessing = time.time()
        print(f"Processing the coalescent tree: {tree_filepath}")
        coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_filepath)
        coalescent_tree.remove_unary_nodes()
        pedigree_filepath = parent_directory_path / pedigree_filename
        print(f"Processing the pedigree: {pedigree_filepath}")
        pedigree_probands = get_pedigree_simulation_probands_for_alignment_mode(coalescent_tree=coalescent_tree)
        pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_filepath,
                                                                             probands=pedigree_probands,
                                                                             preprocess_graph=True)
        end_preprocessing = time.time()
        print(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")
        initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
        alignment_result_path = result_directory_path / tree_filename
        graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                     initial_mapping=initial_mapping, logs_path=alignment_result_path)
        run_alignments_and_save_results(graph_matcher=graph_matcher, output_path=alignment_result_path)


def run_alignments_and_save_results(graph_matcher: GraphMatcher, output_path: str | Path):
    start_alignment = time.time()
    print("Running the alignment")
    result = graph_matcher.find_mapping()
    end_alignment = time.time()
    print(f"Matching time: {end_alignment - start_alignment} seconds")
    graph_matcher.log(f"Matching time: {end_alignment - start_alignment} seconds")
    save_alignment_result_to_files(alignment_result=result, coalescent_tree=graph_matcher.coalescent_tree,
                                   pedigree=graph_matcher.pedigree, directory_path=output_path)


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
        extension = os.path.splitext(filename)[1]
        if extension in extension_skip_list or filename == "SKIP":
            continue
        print(f"{filename}")
        alignment_result_path = result_directory_path / filename
        absolute_path = os.path.abspath(file)
        coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=absolute_path)
        coalescent_tree.remove_unary_nodes()
        initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
        graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                     logs_path=alignment_result_path, initial_mapping=initial_mapping)
        start_alignment = time.time()
        result = graph_matcher.find_mapping()
        end_alignment = time.time()
        print(f"Matching time: {end_alignment - start_alignment} seconds")
        graph_matcher.log(f"Matching time: {end_alignment - start_alignment} seconds")
        total_alignment_time += end_alignment - start_alignment
        count += 1
        save_alignment_result_to_files(alignment_result=result, coalescent_tree=coalescent_tree,
                                       pedigree=pedigree, directory_path=str(alignment_result_path))
    print(f"Average inference time {total_alignment_time / count} seconds")
