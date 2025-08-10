from typing import cast

from alignment.alignment_result import TreeAlignmentResults, SuccessCladeAlignmentResults, AlignmentResult, \
    FailedClimbingCladeAlignmentResults
from alignment.graph_matcher import *
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from scripts.alignment_statistics.calculate_statistics import SuccessCladeAlignmentMetadata, \
    CladeAlignmentStatisticsMetadata, ClimbingFailedCladeAlignmentMetadata
from scripts.utility.basic_utility import *

print_enabled = False
extension_skip_list = [".svg", pedigree_extension]


def save_alignment_result_to_files(alignment_result: TreeAlignmentResults, coalescent_tree: CoalescentTree,
                                   pedigree: PotentialMrcaProcessedGraph, directory_path: str | Path = ""):
    result_parent_directory = Path(directory_path)
    for clade_root, clade_alignment_result in alignment_result.clade_root_to_clade_results.items():
        result_directory_name = f"{clade_root}"
        clade_result_directory_path = result_parent_directory / result_directory_name
        clade_metadata = None
        # Determine the clade alignment result type
        match clade_alignment_result:
            case FailedClimbingCladeAlignmentResults():
                clade_metadata = ClimbingFailedCladeAlignmentMetadata(
                    clade_alignment_result=clade_alignment_result,
                    results_filepath=clade_result_directory_path,
                    clade_root=clade_root,
                    coalescent_tree=coalescent_tree
                )
            case SuccessCladeAlignmentResults():
                clade_alignments_metadata = CladeAlignmentStatisticsMetadata(
                    calculate_similarity=calculate_similarity,
                    calculate_distances_histogram=calculate_distances_histogram,
                    calculate_alignments_likelihoods=calculate_likelihood
                )
                clade_metadata = SuccessCladeAlignmentMetadata(
                    coalescent_tree=coalescent_tree, pedigree=pedigree, clade_root=clade_root,
                    results_filepath=clade_result_directory_path, clade_alignments_metadata=clade_alignments_metadata,
                    clade_alignment_result=clade_alignment_result,
                )
        clade_metadata.save()


def save_alignment_result_and_store_vertex_alignment(
        alignment_result: AlignmentResult,
        general_result: TreeAlignmentResults,
        coalescent_tree: CoalescentTree, directory_path: str | Path):
    assert alignment_result.clade_root is not None
    if isinstance(alignment_result, FailedClimbingAlignmentResult):
        general_result.clade_root_to_clade_results[alignment_result.clade_root] = (
            FailedClimbingCladeAlignmentResults(
                failed_climbing_alignment_info=alignment_result
            )
        )
        return
    if isinstance(alignment_result, FullAlignmentResult):
        # Cast to avoid warnings from IDE
        alignment_result = cast(FullAlignmentResult, alignment_result)
        assert alignment_result.vertex_alignment

        def update_pedigree_vertex_appearance_in_edge_alignment(edge_alignment: dict):
            for edge, path in edge_alignment.items():
                # Process every vertex except the last one (to avoid counting the same vertex twice)
                for vertex in path[:-1]:
                    clade_results.pedigree_vertex_to_edge_alignment_appearance_number[vertex] += 1
            # Count the root separately
            clade_root_candidate = alignment_result[alignment_result.clade_root]
            clade_results.pedigree_vertex_to_edge_alignment_appearance_number[clade_root_candidate] += 1

        # If this is a first alignment for a clade, create the object for collecting the results for the clade
        clade_results: SuccessCladeAlignmentResults = general_result.clade_root_to_clade_results.setdefault(
            alignment_result.clade_root,
            SuccessCladeAlignmentResults(clade_root=alignment_result.clade_root)
        )
        current_clade_alignments = clade_results.alignments
        alignment_index = len(current_clade_alignments)
        current_clade_alignments.append(alignment_result)
        directory_path = Path(directory_path)
        result_dir_path = directory_path / str(alignment_result.clade_root)
        os.makedirs(result_dir_path, exist_ok=True)
        result_filepath = result_dir_path / f"alignment_{alignment_index}"
        alignment_result.save_to_file(filepath=result_filepath, tree=coalescent_tree)
        # Calculate how often a pedigree vertex appears in all the edge alignments if they are present
        if alignment_result.edge_alignments:
            for edge_alignment in alignment_result.edge_alignments:
                update_pedigree_vertex_appearance_in_edge_alignment(edge_alignment)
            clade_results.edge_alignment_total_number += len(alignment_result.edge_alignments)
        # The edge alignments can be heavy, so we can get rid of them after writing the results to the file
        alignment_result.edge_alignments = None
        alignment_result.example_edge_alignment = None


def get_store_and_save_vertex_alignment_callback(coalescent_tree: CoalescentTree, directory_path: str | Path) \
        -> TreeAlignmentResults:
    root_vertices = coalescent_tree.get_founders()
    clade_root_to_clade_results = {
        x: SuccessCladeAlignmentResults(
            clade_root=x
        )
        for x in root_vertices
    }
    alignment_result = TreeAlignmentResults(
        clade_root_to_clade_results=clade_root_to_clade_results
    )
    return alignment_result, partial(save_alignment_result_and_store_vertex_alignment,
                                     coalescent_tree=coalescent_tree,
                                     directory_path=directory_path,
                                     general_result=alignment_result
                                     )


def store_vertex_alignment_callback(alignment_result: AlignmentResult, result: TreeAlignmentResults):
    assert alignment_result.clade_root is not None
    clade_root = alignment_result.clade_root
    match alignment_result:
        case FailedClimbingAlignmentResult():
            result.clade_root_to_clade_results[clade_root] = FailedClimbingCladeAlignmentResults(
                failed_climbing_alignment_info=alignment_result
            )
        case FullAlignmentResult(vertex_alignment=vertex_alignment):
            assert vertex_alignment
            clade_results: SuccessCladeAlignmentResults = result.clade_root_to_clade_results.setdefault(
                clade_root,
                SuccessCladeAlignmentResults(clade_root=clade_root)
            )
            clade_results.alignments.append(alignment_result)


def get_store_vertex_alignment_callback():
    result = TreeAlignmentResults(
        clade_root_to_clade_results={}
    )
    return result, partial(store_vertex_alignment_callback, result_dict=result)


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
        tree_probands = coalescent_tree.get_sink_vertices()
        pedigree_probands = get_pedigree_simulation_probands_for_alignment_mode(vertices=tree_probands)
        pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_filepath,
                                                                             probands=pedigree_probands,
                                                                             preprocess_graph=True)
        end_preprocessing = time.time()
        print(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")
        initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
        alignment_result_path = result_directory_path / tree_filename
        alignment_general_results, save_results_callback = get_store_and_save_vertex_alignment_callback(
            coalescent_tree=coalescent_tree,
            directory_path=alignment_result_path
        )
        graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                     initial_mapping=initial_mapping, logs_path=alignment_result_path,
                                     result_callback_function=save_results_callback)
        run_alignments(graph_matcher=graph_matcher)
        save_alignment_result_to_files(
            alignment_result=alignment_general_results,
            coalescent_tree=coalescent_tree,
            pedigree=pedigree,
            directory_path=alignment_result_path
        )


def run_alignments(graph_matcher: GraphMatcher):
    start_alignment = time.time()
    print("Running the alignment")
    graph_matcher.find_alignments()
    end_alignment = time.time()
    matching_time_message = f"Matching time: {end_alignment - start_alignment} seconds"
    print(matching_time_message)
    graph_matcher.log(matching_time_message)


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
        alignment_general_results, save_results_callback = get_store_and_save_vertex_alignment_callback(
            coalescent_tree=coalescent_tree,
            directory_path=alignment_result_path
        )
        graph_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=pedigree,
                                     logs_path=alignment_result_path, initial_mapping=initial_mapping,
                                     result_callback_function=save_results_callback)
        start_alignment = time.time()
        graph_matcher.find_alignments()
        end_alignment = time.time()
        graph_matcher.log(f"Matching time: {end_alignment - start_alignment} seconds")
        total_alignment_time += end_alignment - start_alignment
        count += 1
        save_alignment_result_to_files(alignment_result=alignment_general_results, coalescent_tree=coalescent_tree,
                                       pedigree=pedigree, directory_path=str(alignment_result_path))
    print(f"Average inference time {total_alignment_time / count} seconds")
