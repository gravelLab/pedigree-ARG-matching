from alignment.graph_matcher import *
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.alignment_statistics.calculate_statistics import CladeMetadata, CladeAlignmentsMetadata
from scripts.utility.basic_utility import *

print_enabled = False
extension_skip_list = [".svg", pedigree_extension]


def save_alignment_result_to_files(alignment_result: {int: [dict]}, coalescent_tree: CoalescentTree,
                                   pedigree: PotentialMrcaProcessedGraph, directory_path: str | Path = ""):
    result_parent_directory = Path(directory_path)
    for clade_root in alignment_result:
        result_directory_name = f"{clade_root}"
        clade_result_directory_path = result_parent_directory / result_directory_name
        valid_alignments = alignment_result[clade_root]
        clade_alignments_metadata = CladeAlignmentsMetadata(
            clade_alignments=valid_alignments,
            calculate_similarity=calculate_similarity,
            calculate_distances_histogram=calculate_distances_histogram,
            calculate_alignments_likelihoods=calculate_likelihood
        )
        clade_metadata = CladeMetadata.get_clade_basic_metadata(
            coalescent_tree=coalescent_tree, pedigree=pedigree, clade_root=clade_root,
            results_filepath=clade_result_directory_path, clade_alignments_metadata=clade_alignments_metadata
        )
        clade_metadata.save()


def process_pedigree_tree_directory(directory: str, result_directory_name: str):
    parent_directory_path = Path(directory)
    result_directory_path = parent_directory_path / result_directory_name
    os.makedirs(result_directory_path)
    directory_files = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]
    pedigree_filename = get_unique_filename_with_specified_extension(directory_path=directory,
                                                                     extension=pedigree_extension)
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
    matching_time_message = f"Matching time: {end_alignment - start_alignment} seconds"
    print(matching_time_message)
    graph_matcher.log(matching_time_message)
    save_alignment_result_to_files(alignment_result=result, coalescent_tree=graph_matcher.coalescent_tree,
                                   pedigree=graph_matcher.pedigree, directory_path=output_path)


def run_alignment_with_multiple_clades_and_save_results(directory: str, result_directory_name: str,
                                                        pedigree: PotentialMrcaProcessedGraph = None):
    parent_directory_path = Path(directory)
    result_directory_path = parent_directory_path / result_directory_name
    os.makedirs(result_directory_path)
    if pedigree is None:
        pedigree_filename = get_unique_filename_with_specified_extension(directory_path=directory,
                                                                         extension=pedigree_extension)
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
