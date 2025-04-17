from __future__ import annotations

import argparse
import csv
import os.path
import sys
import traceback
import warnings
from abc import ABC, abstractmethod
from concurrent.futures import as_completed, ProcessPoolExecutor

from alignment.graph_matcher import *
from scripts import utility
from scripts.run_alignment import save_alignment_result_to_files
from scripts.utility import *

log_directory = "log_dir"
initial_alignment_dir_name = "initial_alignment"
general_result_filename = "simulation_result.txt"
detailed_result_filename = "detailed_simulation_result.txt"
results_csv_filename = "results_by_clade.csv"
partial_results_csv_filename = "partial_results.csv"
csv_header = ("proband_number,no_solutions,individual_and_spouses,individual_and_non_spouse,"
              "no_individual_spouse,neither_individual_nor_spouse,neither_individual_nor_spouse_only_super_founders\n")
total_results_filename = "total_results.txt"
line_content_divider = "------------------------------------\n"
line_section_divider = "#####################################\n"

script_help_description = """
                          This script performs the alignment algorithm between the specified pedigree(s) and
                          coalescent tree(s) and calculates the resulting statistics about these alignments.
                          The primary goal of this script it to automate the alignment process between multiple
                          pedigrees and coalescent trees that have been modified to simulate errors that occur within
                          real data, but it can be used for other purposes as well.

                          You can specify the running mode for the script and the corresponding
                          parameters as the arguments to the script or by running the script
                          without any arguments to start the interactive session.
                          
                          The initial mapping is deducted based on the current configuration in all the modes.
                          
                          Available modes:
                          1) Specify the paths to one pedigree and one coalescent tree on which the alignments should
                          be performed. Then, specify the path to a directory containing directories with pedigree
                          files. You can get this kind of output from the error_pedigree_simulation script.
                          This mode aligns the specified coalescent tree ith all the pedigrees specified
                          in the second option and compares them against the original alignment (the alignment between
                          the tree and the pedigree within the same folder from the first option).
                          
                          2) Specify multiple paths to directories containing a pedigree and a coalescent tree (the 
                          pedigree file must have the .pedigree extension, so that the script can distinguish 
                          between them). Then, specify the path to a directory containing directories with pedigree
                          files.
                          This mode acts as a wrapper for the first mode by doing the same actions for all the trees.
                          
                          3)  Specify the path to a directory containing directories with coalescent trees 
                          and pedigrees.
                          Then, the specify the path to a directory containing directories with pedigree
                          files. This option acts as a wrapper for the second option allowing you to avoid specifying
                          all the paths separately if they directories are within the same folder.
                          """

super_founders = []


def prompt_top_founders():
    global super_founders
    should_process_top_founders = get_yes_or_no("Do you want to specify the top founders?")
    if not should_process_top_founders:
        return
    founders_csv_filepath = get_filepath("Specify the path to the csv file with the top founders:")
    new_super_founders = []
    with open(founders_csv_filepath, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            vertex_id, _ = row[0], float(row[1])
            try:
                vertex_id = int(vertex_id)
            except ValueError:
                warnings.warn(f"The super founder id {vertex_id} is skipped as it's not an integer")
            new_super_founders.append(int(vertex_id))
    super_founders = frozenset(new_super_founders)
    print(f"Identified {len(super_founders)} super founders")


class AlignmentClassification:

    def __init__(self):
        self.total_simulation_number = 0
        self.simulations_with_solutions = 0
        # Keeping track of the number of simulations with 0 solutions
        self.no_solutions_counter = 0
        # Keeping track of the number of simulations where the root assignments are the root vertex
        # itself and their spouses
        # In other words, it is a subset of the correct assignments
        self.only_individual_and_spouses_counter = 0
        # The number of simulations where the correct individual is found, but a non-spouse assignment
        # for the root is also present
        self.individual_and_non_spouse_present_counter = 0
        # The number of simulations where the individual is not present, but one of the spouses is present
        self.individual_not_present_spouse_present_counter = 0
        # The number of simulations where neither the individual nor a spouse is present
        self.neither_individual_nor_spouse_present_counter = 0
        # A subset of the previous case when there is a superfounder assignment
        self.superfounder_present_individual_and_spouses_not_present = 0
        # The next fields keep track of the initial and new assignments to the root vertex
        self.total_number_of_root_vertex_individual_assignments = 0
        self.initial_root_vertex_individual_assignments_number = 0
        self.root_assignment_occurrences = defaultdict(int)

    def _add_root_assigment_dictionary(self, other_root_assignment_occurrences):
        for key, value in other_root_assignment_occurrences.items():
            if key in self.root_assignment_occurrences:
                self.root_assignment_occurrences[key] += value
            else:
                self.root_assignment_occurrences[key] = value

    def add_root_assignments(self, root_assignments: [int]):
        for root_assignment in root_assignments:
            self.root_assignment_occurrences[root_assignment] += 1

    def add_results(self, other_result: AlignmentClassification):
        self.simulations_with_solutions += other_result.simulations_with_solutions
        self.no_solutions_counter += other_result.no_solutions_counter
        self.only_individual_and_spouses_counter += other_result.only_individual_and_spouses_counter
        self.individual_and_non_spouse_present_counter += other_result.individual_and_non_spouse_present_counter
        self.individual_not_present_spouse_present_counter += (
            other_result.individual_not_present_spouse_present_counter)
        self.neither_individual_nor_spouse_present_counter += (
            other_result.neither_individual_nor_spouse_present_counter)
        self.superfounder_present_individual_and_spouses_not_present += (
            other_result.superfounder_present_individual_and_spouses_not_present)
        self.total_number_of_root_vertex_individual_assignments += (
            other_result.total_number_of_root_vertex_individual_assignments)
        self.total_simulation_number += other_result.total_simulation_number
        self._add_root_assigment_dictionary(other_result.root_assignment_occurrences)
        return self

    def csv_row(self, proband_number: str):
        row = (f"{proband_number},{self.no_solutions_counter},"
               f"{self.only_individual_and_spouses_counter},"
               f"{self.individual_and_non_spouse_present_counter},"
               f"{self.individual_not_present_spouse_present_counter},"
               f"{self.neither_individual_nor_spouse_present_counter},"
               f"{self.superfounder_present_individual_and_spouses_not_present}\n")
        return row

    def is_empty(self):
        return self.no_solutions_counter == 0 and self.simulations_with_solutions == 0


def save_multiple_alignment_classifications(results_csv_file,
                                            alignment_classifications: [AlignmentClassification],
                                            proband_number: str):
    for proband_alignment_classification in alignment_classifications:
        row = proband_alignment_classification.csv_row(proband_number=proband_number)
        results_csv_file.write(row)


def save_alignment_results_sorted_by_proband_number(
                        results_csv_filepath: str | Path,
                        proband_number_to_result: dict[int, AlignmentClassification]) -> AlignmentClassification:
    global_alignment_classification = AlignmentClassification()
    with open(results_csv_filepath, "a") as results_csv_file:
        results_csv_file.write(csv_header)
        for proband_number, proband_number_classifications in sorted(proband_number_to_result.items()):
            for proband_number_classification in proband_number_classifications:
                row = proband_number_classification.csv_row(proband_number=str(proband_number))
                results_csv_file.write(row)
                global_alignment_classification.add_results(proband_number_classification)
    return global_alignment_classification


class RootVertexSpouses:
    # This is a simple data class that contains the information about a pedigree ploid vertex,
    # its corresponding individual id, and the spouses

    def __init__(self, vertex_ploid_id: int, individual_and_spouses: set[int] | frozenset[int]):
        self.vertex_ploid_id = vertex_ploid_id
        self.vertex_individual_id = vertex_ploid_id // 2
        if self.vertex_individual_id not in individual_and_spouses:
            raise ValueError(f"The vertex individual id {self.vertex_individual_id} does not belong to "
                             f" the passed set {individual_and_spouses}")
        self.individual_and_spouses = individual_and_spouses

    @staticmethod
    def get_root_vertex_information_from_pedigree(root_vertex_ploid_id: int, pedigree: SimpleGraph):
        root_vertex_children_all_ploids = [
            child_ploid
            for child in pedigree.children_map[root_vertex_ploid_id]
            for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
        ]
        root_vertex_pedigree_individual_and_spouses = {p // 2 for child_ploid
                                                       in root_vertex_children_all_ploids
                                                       for p in pedigree.parents_map.get(child_ploid, [])}
        return RootVertexSpouses(
            vertex_ploid_id=root_vertex_ploid_id,
            individual_and_spouses=root_vertex_pedigree_individual_and_spouses
        )


class AlignmentResultsFileAccess:
    def __init__(self, simulation_directory_path: Path,
                 simulation_general_result_filename: str = general_result_filename,
                 simulation_detailed_result_filename: str = detailed_result_filename,
                 save_alignments: bool = True):
        self.general_result_filepath = simulation_directory_path / simulation_general_result_filename
        self.detailed_result_filepath = simulation_directory_path / simulation_detailed_result_filename
        self.detailed_result_file = None
        self.general_result_file = None
        self.save_alignments = save_alignments

    def open_files(self):
        self.general_result_file = open(self.general_result_filepath, 'a')
        self.detailed_result_file = open(self.detailed_result_filepath, 'a')

    def log_solutions_to_detailed_results_file(self, alignments: list[dict]):
        # Saving the information about the alignments into the detailed file
        self.detailed_result_file.write(line_section_divider)
        self.detailed_result_file.write("Next simulation\n")
        self.detailed_result_file.write(f"Number of alignments in this simulation: {len(alignments)}\n")
        for alignment in alignments:
            self.detailed_result_file.write(line_content_divider)
            self.detailed_result_file.write(f"{alignment}\n")
            self.detailed_result_file.write(line_content_divider)
        self.detailed_result_file.write(line_section_divider)
        self.detailed_result_file.flush()

    def log_solutions_to_general_results_file(self, alignments: list[dict]):
        # Saving the information about the alignments into the general file
        self.general_result_file.write(line_section_divider)
        self.general_result_file.write("Next simulation step\n")
        self.general_result_file.write(
            f"Number of alignments in this simulation {len(alignments)}\n")
        self.general_result_file.flush()

    def log_solutions(self, alignments: list[dict]):
        self.log_solutions_to_detailed_results_file(alignments)
        self.log_solutions_to_general_results_file(alignments)

    def close_files(self):
        if self.general_result_file:
            self.general_result_file.close()
            self.general_result_file = None
        if self.detailed_result_file:
            self.detailed_result_file.close()
            self.detailed_result_file = None


class AbstractAlignmentTask(ABC):

    def __init__(self, file_access: AlignmentResultsFileAccess,
                 alignment_classification: AlignmentClassification = None,
                 root_vertex_info: RootVertexSpouses = None):
        if not alignment_classification:
            alignment_classification = AlignmentClassification()
        self.file_access = file_access
        self.error_results = alignment_classification
        # If root_vertex_info is None, it will be deducted automatically
        self.root_vertex_info = root_vertex_info

    def classify_results(self, coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                         result: dict):
        # The root of the tree can, in principle, be different from the root of the initial, error-free
        # tree. We need to know the id to access the alignments for the clade's root
        coalescent_tree_root = coalescent_tree.get_root_vertex()
        if not self.root_vertex_info:
            # Since the root vertex pedigree candidates are not specified,
            # deduct them automatically (working with msprime simulations)
            self.root_vertex_info = RootVertexSpouses.get_root_vertex_information_from_pedigree(
                pedigree=pedigree,
                root_vertex_ploid_id=coalescent_tree_root
            )
        resulting_alignments = [d for lst in result.values() for d in lst]
        self.error_results.total_simulation_number += 1
        # Get the root assignments for this simulation
        simulation_root_vertex_individual_assignments = {x[coalescent_tree_root] // 2 for x in
                                                         resulting_alignments}
        self.file_access.open_files()
        self.file_access.log_solutions(alignments=resulting_alignments)
        if len(resulting_alignments) > 0:
            self.file_access.general_result_file.write(f"Root vertex individual assignments "
                                                       f"({len(simulation_root_vertex_individual_assignments)}):"
                                                       f"{simulation_root_vertex_individual_assignments}\n")
            # Can close the files earlier to let other processes work with them
            self.file_access.close_files()
            self.error_results.total_number_of_root_vertex_individual_assignments += (
                len(simulation_root_vertex_individual_assignments))
            self.error_results.add_root_assignments(simulation_root_vertex_individual_assignments)
            self.error_results.simulations_with_solutions += 1
            # Classify the simulation result
            if self.root_vertex_info.vertex_individual_id in simulation_root_vertex_individual_assignments:
                if simulation_root_vertex_individual_assignments.issubset(
                        self.root_vertex_info.individual_and_spouses):
                    self.error_results.only_individual_and_spouses_counter += 1
                else:
                    self.error_results.individual_and_non_spouse_present_counter += 1
            else:
                if simulation_root_vertex_individual_assignments.intersection(
                        self.root_vertex_info.individual_and_spouses):
                    self.error_results.individual_not_present_spouse_present_counter += 1
                else:
                    if simulation_root_vertex_individual_assignments.intersection(super_founders):
                        self.error_results.superfounder_present_individual_and_spouses_not_present += 1
                    else:
                        self.error_results.neither_individual_nor_spouse_present_counter += 1
        else:
            self.error_results.no_solutions_counter += 1
        self.file_access.close_files()

    @abstractmethod
    def run(self):
        pass


class AlignmentTask(AbstractAlignmentTask):
    def __init__(self, pedigree_path: Path, coalescent_tree_path: Path,
                 file_access: AlignmentResultsFileAccess,
                 simulation_subdirectory_path: str | Path,
                 root_vertex_info: RootVertexSpouses = None,
                 alignment_classification: AlignmentClassification = None):
        super().__init__(file_access=file_access, alignment_classification=alignment_classification,
                         root_vertex_info=root_vertex_info)
        self.pedigree_path = pedigree_path
        self.coalescent_tree_path = coalescent_tree_path
        self.file_access = file_access
        self.simulation_subdirectory_path = Path(simulation_subdirectory_path)

    def run(self):
        print(f"Parsing {self.pedigree_path} [{self.coalescent_tree_path}]")
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=self.coalescent_tree_path)
        coalescent_tree.remove_unary_nodes()
        ascending_genealogy_probands = get_pedigree_simulation_probands_for_alignment_mode(
            coalescent_tree=coalescent_tree
        )
        pedigree = (PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=self.pedigree_path,
                                                                              preprocess_graph=True,
                                                                              probands=ascending_genealogy_probands
                                                                              )
                    )
        if self.file_access.save_alignments:
            os.mkdir(self.simulation_subdirectory_path)
            log_directory_path = str(self.simulation_subdirectory_path / logs_default_directory_name)
        else:
            log_directory_path = None
        initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
        matcher = GraphMatcher(coalescent_tree=coalescent_tree,
                               processed_graph=pedigree,
                               logs_path=log_directory_path,
                               initial_mapping=initial_mapping)
        result = matcher.find_mapping()
        if self.file_access.save_alignments:
            alignments_dir_path = str(self.simulation_subdirectory_path)
            save_alignment_result_to_files(alignment_result=result,
                                           coalescent_tree=coalescent_tree,
                                           pedigree=pedigree,
                                           directory_path=alignments_dir_path)
        self.classify_results(result=result, coalescent_tree=coalescent_tree,
                              pedigree=pedigree)
        return self.error_results


class MaximumAlignableSubcladeTask(AbstractAlignmentTask):

    def __init__(self, pedigree: PotentialMrcaProcessedGraph, tree_cut_clades_dir_path: str | Path,
                 file_access: AlignmentResultsFileAccess,
                 simulation_subdirectory_path: str | Path,
                 root_vertex_info: RootVertexSpouses = None,
                 alignment_classification: AlignmentClassification = None
                 ):
        super().__init__(file_access=file_access, alignment_classification=alignment_classification,
                         root_vertex_info=root_vertex_info)
        self.pedigree = pedigree
        self.tree_subclades_dir_path = Path(tree_cut_clades_dir_path)
        self.simulation_subdirectory_path = Path(simulation_subdirectory_path)

    def process_tree_subclades_directory(self):
        tree_files = [file for file in self.tree_subclades_dir_path.iterdir() if file.is_file()]
        # Extract filenames and check if they are valid natural numbers
        file_numbers = []
        for file in tree_files:
            if file.name.startswith('.'):
                continue
            try:
                num = int(file.stem)
                if num < 0:
                    raise ValueError(f"Invalid filename: {file.name}. Expected non-negative integers.")
                file_numbers.append((num, file))
            except ValueError:
                raise ValueError(f"Invalid filename: {file.name}. Expected numeric filenames.")
        # Sort by the numeric value of the filename
        file_numbers.sort(key=lambda x: x[0])
        # Ensure all numbers from 0 to max exist
        expected_numbers = set(range(len(file_numbers)))
        actual_numbers = {num for num, _ in file_numbers}
        if expected_numbers != actual_numbers:
            return None
        return file_numbers

    def run(self) -> (int, AlignmentClassification):
        proband_number = None
        processed_tree_subclade_paths = self.process_tree_subclades_directory()
        if not processed_tree_subclade_paths:
            warnings.warn(f"Skipping, received an invalid tree-subclades directory path: "
                          f"{self.tree_subclades_dir_path}")
            return proband_number, self.error_results
        for tree_index, tree_path in processed_tree_subclade_paths:
            coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path)
            initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
            tree_log_directory = (self.simulation_subdirectory_path / f"{tree_index}"
                                  / logs_default_directory_name)
            os.makedirs(tree_log_directory, exist_ok=True)
            matcher = GraphMatcher(coalescent_tree=coalescent_tree,
                                   processed_graph=self.pedigree,
                                   logs_path=tree_log_directory,
                                   initial_mapping=initial_mapping)
            result = matcher.find_mapping()
            resulting_alignments = [d for lst in result.values() for d in lst]
            if resulting_alignments:
                # If the clade is alignable, stop
                proband_number = len(coalescent_tree.get_probands())
                self.classify_results(coalescent_tree=coalescent_tree, pedigree=self.pedigree,
                                      result=result)
                break
        return proband_number, self.error_results


def run_initial_alignment(coalescent_tree: CoalescentTree, initial_pedigree: PotentialMrcaProcessedGraph,
                          result_filepath: str, save_alignments_to_files: bool):
    result_filepath = Path(result_filepath)
    os.mkdir(result_filepath)
    if save_alignments_to_files:
        log_directory_path = result_filepath / logs_default_directory_name
    else:
        log_directory_path = None
    initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
    initial_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=initial_pedigree,
                                   initial_mapping=initial_mapping, logs_path=log_directory_path)
    initial_result_dict = initial_matcher.find_mapping()
    if save_alignments_to_files:
        save_alignment_result_to_files(alignment_result=initial_result_dict,
                                       coalescent_tree=coalescent_tree,
                                       pedigree=initial_pedigree,
                                       directory_path=result_filepath)
    return initial_result_dict


class BaseErrorAlignmentComparison(ABC):

    def __init__(self, initial_coalescent_tree_path: str, initial_pedigree_path: str,
                 simulation_folder_path: str | Path,
                 save_alignments_to_files: bool = True):
        self.error_comparison_results = AlignmentClassification()
        self.initial_coalescent_tree_path = initial_coalescent_tree_path
        self.initial_pedigree_path = initial_pedigree_path
        self.save_alignments_to_files = save_alignments_to_files
        self.simulation_folder_directory = Path(simulation_folder_path)
        self.simulation_folder_directory.mkdir(parents=True, exist_ok=True)
        self.file_access = AlignmentResultsFileAccess(
            simulation_directory_path=self.simulation_folder_directory
        )
        self.log_graph_paths()

    def log_graph_paths(self):
        general_result_file = open(self.file_access.general_result_filepath, 'a')
        general_result_file.write(f"The path to the initial, error-free pedigree is:"
                                  f" {self.initial_pedigree_path}\n")
        general_result_file.write(f"The path to the coalescent tree is: {self.initial_coalescent_tree_path}\n")
        general_result_file.write(line_section_divider)
        general_result_file.close()

    def log_overall_results(self):
        self.error_comparison_results.simulations_with_solutions = (
                self.error_comparison_results.total_simulation_number -
                self.error_comparison_results.no_solutions_counter)
        general_result_file = open(self.file_access.general_result_filepath, 'a')
        general_result_file.write(line_content_divider)
        general_result_file.write(f"The number of simulations with solutions: "
                                  f"{self.error_comparison_results.simulations_with_solutions}\n")
        general_result_file.write(f"The number of simulations with no solutions: "
                                  f"{self.error_comparison_results.no_solutions_counter}\n")
        general_result_file.write(f"The number of simulations with only the individual and spouses: "
                                  f"{self.error_comparison_results.only_individual_and_spouses_counter}\n")
        general_result_file.write(f"The number of simulations with the individual and a non-spouse present: "
                                  f"{self.error_comparison_results.individual_and_non_spouse_present_counter}\n")
        general_result_file.write(f"The number of simulations without the individual, but with a spouse: "
                                  f"{self.error_comparison_results.individual_not_present_spouse_present_counter}\n")
        general_result_file.write(f"The number of simulations without the individual and all the spouses "
                                  f"(at least 1 non-superfounder): "
                                  f"{self.error_comparison_results.neither_individual_nor_spouse_present_counter}\n")
        general_result_file.write(
            f"The number of simulations without the individual and all the spouses "
            f"where there is a superfounder: "
            f"{self.error_comparison_results.superfounder_present_individual_and_spouses_not_present}\n")
        general_result_file.write(line_content_divider)
        general_result_file.write(
            f"Initial number of individual candidates for the root: "
            f"{self.error_comparison_results.initial_root_vertex_individual_assignments_number}\n"
        )
        if self.error_comparison_results.simulations_with_solutions != 0:
            average_number_candidates = (
                    self.error_comparison_results.total_number_of_root_vertex_individual_assignments
                    / self.error_comparison_results.simulations_with_solutions)
            general_result_file.write(
                f"Average number of individual candidates for the root per simulation with solutions: "
                f"{average_number_candidates}\n")
        general_result_file.close()

    def run_initial_alignment_and_log_general_graph_information(self) -> RootVertexSpouses:
        initial_alignment_result_path = self.simulation_folder_directory / initial_alignment_dir_name
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(
            filepath=self.initial_coalescent_tree_path
        )
        coalescent_tree.remove_unary_nodes()
        number_of_probands = len(coalescent_tree.probands)
        non_proband_vertices = len(coalescent_tree.vertex_to_level_map) - number_of_probands
        general_result_file = open(self.file_access.general_result_filepath, 'a')
        general_result_file.write(f"The number of probands is {number_of_probands}\n")
        general_result_file.write(f"The number of non-proband vertices is: {non_proband_vertices}\n")
        general_result_file.write(line_content_divider)
        ascending_genealogy_probands = get_pedigree_simulation_probands_for_alignment_mode(
            coalescent_tree=coalescent_tree
        )
        should_run_initial_alignment = False
        initial_pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(
            filepath=self.initial_pedigree_path,
            probands=ascending_genealogy_probands,
            preprocess_graph=should_run_initial_alignment
        )
        root_vertex = coalescent_tree.get_root_vertex()
        root_vertex_info = RootVertexSpouses.get_root_vertex_information_from_pedigree(
            root_vertex_ploid_id=root_vertex, pedigree=initial_pedigree
        )
        assert root_vertex_info.individual_and_spouses, "Invalid pedigree"
        if should_run_initial_alignment:
            initial_result_dict = run_initial_alignment(coalescent_tree=coalescent_tree,
                                                        initial_pedigree=initial_pedigree,
                                                        save_alignments_to_files=self.save_alignments_to_files,
                                                        result_filepath=initial_alignment_result_path)
            valid_alignments = [d for lst in initial_result_dict.values() for d in lst]
            # Log the general information about the graphs and the initial alignment
            general_result_file.write(f"Total number of initial alignments is: {len(valid_alignments)}\n")
            # Calculating the number of individuals to which the clade's root is assigned
            initial_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in valid_alignments}
            self.error_comparison_results.initial_root_vertex_individual_assignments_number = (
                len(initial_root_vertex_individual_assignments)
            )
            general_result_file.write(
                f"Number of individual assignments for the root: "
                f"{self.error_comparison_results.initial_root_vertex_individual_assignments_number}\n"
            )
            general_result_file.write(
                f"Individual assignments for the root: {initial_root_vertex_individual_assignments}\n")
        general_result_file.write(f"Clade root individual: {root_vertex_info.vertex_individual_id}\n")
        general_result_file.write(f"Clade root and their spouses: {root_vertex_info.individual_and_spouses}\n")
        general_result_file.flush()
        general_result_file.close()
        return root_vertex_info

    def run_alignments(self) -> AlignmentClassification:
        alignment_tasks: [AlignmentTask] = self.get_alignment_tasks()
        try:
            for alignment_task in alignment_tasks:
                alignment_task.run()
        except KeyboardInterrupt:
            print("Stop the simulation, log the final results")
        except Exception as ex:
            print(f"Exception occurred: {ex}")
            traceback.print_exc()
        self.log_overall_results()
        return self.error_comparison_results

    @abstractmethod
    def get_alignment_tasks(self):
        pass


class ErrorPedigreesAlignmentComparison(BaseErrorAlignmentComparison):

    def __init__(self, coalescent_tree_path: str, error_free_pedigree_path: str,
                 error_pedigrees_folder: str, simulation_folder_subpath: str,
                 save_alignments_to_files: bool = True):
        error_pedigree_parent_directory = Path(error_pedigrees_folder).parent
        simulation_directory_path = error_pedigree_parent_directory / simulation_folder_subpath
        super().__init__(
            initial_coalescent_tree_path=coalescent_tree_path,
            initial_pedigree_path=error_free_pedigree_path,
            simulation_folder_path=simulation_directory_path,
            save_alignments_to_files=save_alignments_to_files
        )
        self.error_pedigrees_folder = error_pedigrees_folder

    def get_alignment_tasks(self) -> [AlignmentTask]:
        root_vertex_info = self.run_initial_alignment_and_log_general_graph_information()
        current_directory = Path(self.error_pedigrees_folder)
        alignment_tasks = []
        for subdirectory in os.listdir(self.error_pedigrees_folder):
            simulation_subdirectory_path = self.simulation_folder_directory / subdirectory
            subdirectory_path = current_directory / subdirectory
            pedigree_file = get_unique_filename_with_specified_extension(directory_path=subdirectory_path,
                                                                         extension=".pedigree")

            pedigree_filepath = subdirectory_path / pedigree_file
            alignment_task = AlignmentTask(
                pedigree_path=pedigree_filepath,
                coalescent_tree_path=self.initial_coalescent_tree_path,
                file_access=self.file_access,
                simulation_subdirectory_path=simulation_subdirectory_path,
                alignment_classification=self.error_comparison_results,
                root_vertex_info=root_vertex_info
            )
            alignment_tasks.append(alignment_task)
        return alignment_tasks


class ErrorTreesAlignmentComparison(BaseErrorAlignmentComparison):

    def __init__(self, initial_tree_path: str, initial_pedigree_path: str,
                 tree_error_directory_path: str,
                 simulation_folder_subpath: str, save_alignments_to_files: bool = True):
        error_pedigree_parent_directory = Path(tree_error_directory_path).parent
        simulation_directory_path = error_pedigree_parent_directory / simulation_folder_subpath
        super().__init__(
            initial_pedigree_path=initial_pedigree_path,
            initial_coalescent_tree_path=initial_tree_path,
            save_alignments_to_files=save_alignments_to_files,
            simulation_folder_path=simulation_directory_path
        )
        self.tree_error_directory_path = tree_error_directory_path

    def get_alignment_tasks(self) -> [AlignmentTask]:
        root_vertex_info = self.run_initial_alignment_and_log_general_graph_information()
        # Loop over the error trees and generate the alignment tasks
        current_directory = Path(self.tree_error_directory_path)
        alignment_tasks = []
        for coalescent_tree_file in os.listdir(self.tree_error_directory_path):
            coalescent_tree_filepath = current_directory / coalescent_tree_file
            if not os.path.isfile(coalescent_tree_filepath):
                continue
            simulation_subdirectory_path = self.simulation_folder_directory / coalescent_tree_file
            alignment_task = AlignmentTask(
                pedigree_path=self.initial_pedigree_path,
                coalescent_tree_path=coalescent_tree_filepath,
                file_access=self.file_access,
                simulation_subdirectory_path=simulation_subdirectory_path,
                alignment_classification=self.error_comparison_results,
                root_vertex_info=root_vertex_info
            )
            alignment_tasks.append(alignment_task)
        return alignment_tasks


def run_single_data_pedigree_error_directory_alignment(tree_path: str, pedigree_no_errors_path: str,
                                                       error_pedigrees_folder_path: str, input_simulation_subpath: str):
    os.chdir(error_pedigrees_folder_path)
    os.chdir("..")
    error_pedigree_alignment_comparison = ErrorPedigreesAlignmentComparison(
        coalescent_tree_path=tree_path,
        error_pedigrees_folder=error_pedigrees_folder_path,
        error_free_pedigree_path=pedigree_no_errors_path,
        simulation_folder_subpath=input_simulation_subpath)
    return error_pedigree_alignment_comparison.run_alignments()


def run_single_data_tree_error_directory_alignment(tree_path: str, pedigree_path: str,
                                                   error_trees_folder_path: str, simulation_subpath: str) \
        -> AlignmentClassification:
    os.chdir(error_trees_folder_path)
    os.chdir("..")
    error_pedigree_alignment_comparison = ErrorTreesAlignmentComparison(
        initial_tree_path=tree_path,
        initial_pedigree_path=pedigree_path,
        tree_error_directory_path=error_trees_folder_path,
        simulation_folder_subpath=simulation_subpath
    )
    return error_pedigree_alignment_comparison.run_alignments()


def single_initial_data_directory_pedigree_error_directory_alignment_mode():
    current_working_directory = Path.cwd()
    tree_path = get_filepath("Specify the absolute path to the initial coalescent tree:")
    pedigree_no_errors_path = get_filepath("Specify the absolute path to the initial pedigree:")
    error_pedigrees_folder_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    os.chdir(error_pedigrees_folder_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    os.chdir(current_working_directory)
    run_single_data_pedigree_error_directory_alignment(tree_path=tree_path,
                                                       pedigree_no_errors_path=pedigree_no_errors_path,
                                                       error_pedigrees_folder_path=error_pedigrees_folder_path,
                                                       input_simulation_subpath=input_simulation_subpath)


def single_initial_data_directory_tree_error_directory_alignment_mode():
    current_working_directory = Path.cwd()
    tree_path = get_filepath("Specify the absolute path to the initial coalescent tree:")
    pedigree_path = get_filepath("Specify the absolute path to the initial pedigree:")
    error_trees_folder_path = get_directory_path("Specify the absolute path to the error trees directory:")
    os.chdir(error_trees_folder_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    os.chdir(current_working_directory)
    run_single_data_tree_error_directory_alignment(
        tree_path=tree_path,
        pedigree_path=pedigree_path,
        error_trees_folder_path=error_trees_folder_path,
        simulation_subpath=input_simulation_subpath
    )


def get_unique_folder(base_name: str, directory: Path = Path.cwd()) -> Path:
    # Initial attempt with the base name
    unique_name = base_name
    count = 1

    # Check if the directory exists, and add a prefix if it does
    while (directory / unique_name).exists():
        unique_name = f"{count}_{base_name}"
        count += 1
    return unique_name


def process_path(path_to_process: str | Path, simulation_folder_path: str | Path,
                 error_pedigrees_folder_path: str | Path, simulation_subpath: str) \
        -> ErrorPedigreesAlignmentComparison:
    absolute_path_to_process = simulation_folder_path / path_to_process
    paths = get_paths_from_tree_pedigree_directory(absolute_path_to_process)
    if not paths:
        warnings.warn(f"The {path_to_process} directory is invalid. Skipping.")
        return None
    pedigree_path, tree_path = paths
    simulation_folder_basename = os.path.basename(path_to_process)
    simulation_folder_unique_name = get_unique_folder(base_name=simulation_folder_basename,
                                                      directory=simulation_folder_path)
    simulation_path = str(Path(simulation_subpath) / Path(simulation_folder_unique_name))
    error_pedigree_alignment_comparison = ErrorPedigreesAlignmentComparison(
        coalescent_tree_path=tree_path,
        error_pedigrees_folder=error_pedigrees_folder_path,
        error_free_pedigree_path=pedigree_path,
        simulation_folder_subpath=simulation_path
    )
    return error_pedigree_alignment_comparison


def save_error_alignments_to_file(filepath: str | Path, error_alignments_results: AlignmentClassification):
    with open(filepath, "a") as results_file:
        results_file.write(f"The number of simulations with solutions: "
                           f"{error_alignments_results.simulations_with_solutions}\n")
        results_file.write(f"The number of simulations with no solutions: "
                           f"{error_alignments_results.no_solutions_counter}\n")
        results_file.write(f"The number of simulations with only the individual and spouses: "
                           f"{error_alignments_results.only_individual_and_spouses_counter}\n")
        results_file.write(f"The number of simulations with the individual and a non-spouse present: "
                           f"{error_alignments_results.individual_and_non_spouse_present_counter}\n")
        results_file.write(f"The number of simulations without the individual, but with a spouse: "
                           f"{error_alignments_results.individual_not_present_spouse_present_counter}\n")
        results_file.write(f"The number of simulations without the individual and spouses: "
                           f"{error_alignments_results.neither_individual_nor_spouse_present_counter}\n")
        results_file.write(f"The number of simulations without individual and spouses where the rest "
                           f"are super founders: "
                           f"{error_alignments_results.superfounder_present_individual_and_spouses_not_present}\n")


def save_candidate_appearances(filepath: str | Path, error_alignments_results: AlignmentClassification):
    with open(filepath, "a") as results_file:
        results_file.write(line_section_divider)
        results_file.write("The number of root assignment candidates' occurrences\n"
                           "The format is "
                           "'{root_assignment}: {occurrences} ({occurrences / number_simulations_with_solutions})':\n")
        sorted_assignment_occurrences_dict = dict(sorted(error_alignments_results.root_assignment_occurrences.items(),
                                                         key=lambda item: item[1], reverse=True))
        for root_assignment, occurrences in sorted_assignment_occurrences_dict.items():
            if error_alignments_results.simulations_with_solutions == 0:
                occurrences_portion = 0
            else:
                occurrences_portion = occurrences / error_alignments_results.simulations_with_solutions
            results_file.write(f"{root_assignment}: {occurrences} ({occurrences_portion})\n")


def run_specific_error_pedigree_directories_iterative(paths: list[str | Path], error_pedigrees_folder_path: str,
                                                      simulation_subpath: str | Path):
    error_pedigree_parent_directory = Path(error_pedigrees_folder_path).parent
    simulation_folder_path = error_pedigree_parent_directory / simulation_subpath
    os.makedirs(simulation_folder_path, exist_ok=True)
    # The list will contain tuples of form (proband_number, alignment_classification)
    total_results = []
    for path_to_process in paths:
        error_alignment = process_path(path_to_process=path_to_process,
                                       simulation_folder_path=simulation_folder_path,
                                       error_pedigrees_folder_path=error_pedigrees_folder_path,
                                       simulation_subpath=simulation_subpath)
        if not error_alignment:
            continue
        parent_dir_name = Path(path_to_process).parent.name
        tree_results = error_alignment.run_alignments()
        total_results.append((parent_dir_name, tree_results))
    return total_results


def get_alignments_tasks(error_alignment_comparison: ErrorPedigreesAlignmentComparison):
    return error_alignment_comparison, error_alignment_comparison.get_alignment_tasks()


def process_alignment_task(path_alignment_task: AlignmentTask):
    path_alignment_task.run()
    return path_alignment_task


def run_specific_error_pedigree_directories_parallel(paths: list[str], error_pedigrees_folder_path: str,
                                                     simulation_subpath: str | Path, max_workers: int):
    error_pedigrees_folder_path = Path(error_pedigrees_folder_path)
    error_pedigree_parent_directory = error_pedigrees_folder_path.parent
    simulation_directory_path = error_pedigree_parent_directory / simulation_subpath
    os.makedirs(name=simulation_directory_path, exist_ok=True)

    tree_alignment_comparisons: [ErrorPedigreesAlignmentComparison] = []
    # Separating all the pedigree-tree directories into separate comparison objects
    for path in paths:
        tree_alignment_comparison: ErrorPedigreesAlignmentComparison = process_path(path, simulation_directory_path,
                                                                                    error_pedigrees_folder_path,
                                                                                    simulation_subpath)
        if not tree_alignment_comparison:
            continue
        tree_alignment_comparisons.append(tree_alignment_comparison)
    alignment_tasks: [AlignmentTask] = []
    coalescent_path_to_alignment_comparison: {str: ErrorPedigreesAlignmentComparison} = \
        {x.initial_coalescent_tree_path: x for x in tree_alignment_comparisons}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(get_alignments_tasks, tree_alignments_comparison)
            for tree_alignments_comparison in tree_alignment_comparisons]
        for future in as_completed(futures):
            _, path_alignment_tasks = future.result()
            alignment_tasks.extend(path_alignment_tasks)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_alignment_task, alignment_task)
            for alignment_task in alignment_tasks]
        for future in as_completed(futures):
            alignment_task = future.result()
            alignment_comparison = coalescent_path_to_alignment_comparison[alignment_task.coalescent_tree_path]
            alignment_comparison.error_comparison_results.add_results(alignment_task.error_results)
            print(f"Finished alignment {alignment_task.coalescent_tree_path} | {alignment_task.pedigree_path}")
    # Save and return the results
    clade_results_list = []
    for tree_alignment_comparison in tree_alignment_comparisons:
        tree_alignment_comparison.log_overall_results()
        clade_results = tree_alignment_comparison.error_comparison_results
        proband_number = Path(tree_alignment_comparison.initial_coalescent_tree_path).parent.parent.name
        clade_results_list.append((proband_number, clade_results))
    return clade_results_list


def run_specific_error_pedigree_directories(paths: list[str], error_pedigrees_folder_path: str,
                                            simulation_subpath: str, parallelize: bool = True,
                                            max_workers: int = 0):
    if parallelize:
        results = run_specific_error_pedigree_directories_parallel(
            paths=paths,
            error_pedigrees_folder_path=error_pedigrees_folder_path,
            simulation_subpath=simulation_subpath,
            max_workers=max_workers
        )
    else:
        results = run_specific_error_pedigree_directories_iterative(
            paths=paths,
            error_pedigrees_folder_path=error_pedigrees_folder_path,
            simulation_subpath=simulation_subpath
        )
    error_pedigrees_folder_path = Path(error_pedigrees_folder_path)
    simulation_folder_path = error_pedigrees_folder_path.parent / simulation_subpath
    results_csv_filepath = simulation_folder_path / results_csv_filename
    with open(results_csv_filepath, "a") as results_file:
        results_file.write(csv_header)
        for parent_dir_name, tree_results in results:
            results_file.write(tree_results.csv_row(proband_number=parent_dir_name))
    return results


def specific_error_pedigree_directories_mode():
    current_working_directory = Path.cwd()
    paths = get_directory_paths("Specify the absolute paths to all the tree-pedigree directories")
    error_pedigrees_directory_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    os.chdir(error_pedigrees_directory_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    os.chdir(current_working_directory)
    run_specific_error_pedigree_directories(paths=paths,
                                            error_pedigrees_folder_path=error_pedigrees_directory_path,
                                            simulation_subpath=input_simulation_subpath)


def get_absolute_paths_to_subfolders(directory_path: str) -> list[str]:
    return [
        os.path.abspath(os.path.join(directory_path, subdir))
        for subdir in os.listdir(directory_path)
        if os.path.isdir(os.path.join(directory_path, subdir))
    ]


def children_tree_directories_error_pedigree_directory_mode():
    current_working_directory = Path.cwd()
    error_pedigrees_directory_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    tree_pedigree_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory "
                                                             "containing tree-pedigree subdirectories:")
    os.chdir(error_pedigrees_directory_path)
    input_simulation_name = get_non_existing_path("Specify the simulation subpath:")
    tree_pedigree_paths = get_absolute_paths_to_subfolders(tree_pedigree_parent_directory_path)
    os.chdir(current_working_directory)
    parallel_processing = get_yes_or_no("Do you want to parallelize the alignments?")
    max_workers = 0
    if parallel_processing:
        max_workers = get_natural_number_with_lower_bound(input_request="Specify the number of worker"
                                                                        " processed to be used (at least 2): ",
                                                          lower_bound=2
                                                          )
    run_specific_error_pedigree_directories(paths=tree_pedigree_paths,
                                            error_pedigrees_folder_path=error_pedigrees_directory_path,
                                            simulation_subpath=input_simulation_name,
                                            parallelize=parallel_processing,
                                            max_workers=max_workers)


def process_multiple_tree_pedigree_directories_and_corresponding_tree_error_directories(
        tree_pedigree_parent_directory_path: str | Path, error_trees_parent_directory_path: str | Path,
        input_simulation_subpath: str | Path):
    tree_pedigree_parent_directory_path = Path(tree_pedigree_parent_directory_path)
    error_trees_parent_directory_path = Path(error_trees_parent_directory_path)
    input_simulation_subpath = Path(input_simulation_subpath)
    total_error_comparison_results = []
    for tree_pedigree_directory in os.listdir(tree_pedigree_parent_directory_path):
        tree_pedigree_directory_path = tree_pedigree_parent_directory_path / tree_pedigree_directory
        if not os.path.isdir(tree_pedigree_directory_path):
            continue
        paths = get_paths_from_tree_pedigree_directory(tree_pedigree_directory_path)
        if not paths:
            warnings.warn(f"The subdirectory {tree_pedigree_directory} is not a valid tree-pedigree directory."
                          f"Skipping this simulation")
            continue
        pedigree_path, tree_path = paths
        corresponding_error_trees_directory_path = error_trees_parent_directory_path / tree_pedigree_directory
        if not os.path.isdir(corresponding_error_trees_directory_path):
            warnings.warn(f"The directory with modified trees {corresponding_error_trees_directory_path} which is"
                          f"supposed to correspond to the {tree_pedigree_directory_path} tree-pedigree directory"
                          f" was not found. Skipping this simulation")
            continue
        tree_pedigree_simulation_subpath = input_simulation_subpath / tree_pedigree_directory
        tree_pedigree_error_comparison = run_single_data_tree_error_directory_alignment(
            tree_path=tree_path,
            pedigree_path=pedigree_path,
            error_trees_folder_path=corresponding_error_trees_directory_path,
            simulation_subpath=tree_pedigree_simulation_subpath
        )
        if not tree_pedigree_error_comparison.is_empty():
            total_error_comparison_results.append(tree_pedigree_error_comparison)
    return total_error_comparison_results


def multiple_tree_pedigree_directories_and_corresponding_tree_error_directories():
    tree_pedigree_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory "
                                                             "containing tree-pedigree subdirectories:")
    tree_pedigree_parent_directory_path = Path(tree_pedigree_parent_directory_path)
    error_trees_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory containing"
                                                           " error trees subdirectories:")
    error_trees_parent_directory_path = Path(error_trees_parent_directory_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    input_simulation_subpath = Path(input_simulation_subpath)
    simulation_path = error_trees_parent_directory_path / input_simulation_subpath
    os.makedirs(simulation_path, exist_ok=True)
    total_error_comparison_results = (
        process_multiple_tree_pedigree_directories_and_corresponding_tree_error_directories(
            tree_pedigree_parent_directory_path=tree_pedigree_parent_directory_path,
            error_trees_parent_directory_path=error_trees_parent_directory_path,
            input_simulation_subpath=input_simulation_subpath
        ))
    total_results_filepath = simulation_path / total_results_filename
    save_error_alignments_to_file(total_results_filepath, total_error_comparison_results)
    return total_error_comparison_results


def create_or_resolve_polytomy_script():
    proband_number_parent_directory_path = get_directory_path("Specify the path to the parent directory with "
                                                              "error simulated trees:")
    proband_number_parent_directory_path = Path(proband_number_parent_directory_path)
    initial_data_parent_directory_path = get_directory_path("Specify the path to the parent directory with the "
                                                            "initial data (tree-pedigree pairs grouped "
                                                            "by the proband number): ")
    initial_data_parent_directory_path = Path(initial_data_parent_directory_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    simulation_path = proband_number_parent_directory_path / input_simulation_subpath
    os.mkdir(simulation_path)
    partial_csv_filepath = simulation_path / partial_results_csv_filename
    proband_number_to_result = dict()
    with open(partial_csv_filepath, "a") as partial_result_csv_file:
        partial_result_csv_file.write(csv_header)
        for proband_number_directory in os.listdir(initial_data_parent_directory_path):
            try:
                proband_number = int(proband_number_directory)
            except ValueError:
                warnings.warn(f"Skipping {proband_number_directory}")
            initial_data_proband_dir_path = initial_data_parent_directory_path / proband_number_directory
            corresponding_proband_number_error_trees = proband_number_parent_directory_path / proband_number_directory
            # If the file is not a directory, or it's empty, simply skip the file
            if not os.path.isdir(corresponding_proband_number_error_trees) or not os.listdir(
                    corresponding_proband_number_error_trees):
                continue
            proband_number_classifications = (
                process_multiple_tree_pedigree_directories_and_corresponding_tree_error_directories(
                    tree_pedigree_parent_directory_path=initial_data_proband_dir_path,
                    error_trees_parent_directory_path=corresponding_proband_number_error_trees,
                    input_simulation_subpath=input_simulation_subpath,
                )
            )
            proband_number_to_result[proband_number] = proband_number_classifications
            for proband_number_classification in proband_number_classifications:
                row = proband_number_classification.csv_row(proband_number=str(proband_number))
                partial_result_csv_file.write(row)
            partial_result_csv_file.flush()
    results_csv_filepath = simulation_path / results_csv_filename
    global_alignment_classification = save_alignment_results_sorted_by_proband_number(
                                            results_csv_filepath=results_csv_filepath,
                                            proband_number_to_result=proband_number_classification
    )
    os.remove(partial_csv_filepath)
    return global_alignment_classification


def tree_pedigree_subdirectories():
    parent_directory = Path(get_directory_path("Specify the path to the parent directory:"))
    simulation_name = get_non_empty_string("Specify the simulation name:")
    simulation_directory_path = parent_directory / simulation_name
    os.mkdir(simulation_directory_path)
    alignment_tasks = []
    file_access = AlignmentResultsFileAccess(simulation_directory_path=simulation_directory_path)
    proband_number_to_results = dict()
    for proband_number_directory in os.listdir(parent_directory):
        proband_number_directory_path = parent_directory / proband_number_directory
        if not os.path.isdir(proband_number_directory_path):
            continue
        proband_alignment_classifications = []
        for tree_pedigree_subdirectory in os.listdir(proband_number_directory_path):
            tree_pedigree_subdirectory_path = proband_number_directory_path / tree_pedigree_subdirectory
            paths = get_paths_from_tree_pedigree_directory(tree_pedigree_subdirectory_path)
            if not paths:
                continue
            pedigree_path, tree_path = paths
            proband_alignment_classification = AlignmentClassification()
            simulation_subpath = tree_pedigree_subdirectory_path / simulation_name
            alignment_task = AlignmentTask(pedigree_path=pedigree_path,
                                           coalescent_tree_path=tree_path,
                                           file_access=file_access,
                                           alignment_classification=proband_alignment_classification,
                                           simulation_subdirectory_path=simulation_subpath,
                                           # Deduct root vertex information automatically, as it stays the same
                                           root_vertex_info=None)
            alignment_tasks.append(alignment_task)
            proband_alignment_classifications.append(proband_alignment_classification)
        if not proband_alignment_classifications:
            continue
        proband_number_simulation_subpath = proband_number_directory_path / simulation_name
        os.mkdir(proband_number_simulation_subpath)
        proband_number_to_results[proband_number_simulation_subpath] = proband_alignment_classifications
    for alignment_task in alignment_tasks:
        alignment_task.run()
    alignment_classifications = []
    results_csv_filepath = simulation_directory_path / results_csv_filename
    with open(results_csv_filepath, "a") as results_csv_file:
        results_csv_file.write(csv_header)
        for proband_simulation_path, proband_alignment_classifications in proband_number_to_results.items():
            initial_directory_name = proband_simulation_path.parent.name
            proband_results_path = proband_simulation_path / total_results_filename
            # Calculating the accumulative result for this proband number
            proband_result = reduce(lambda x, y: x.add_results(y), proband_alignment_classifications)
            save_error_alignments_to_file(proband_results_path, proband_result)
            alignment_classifications.extend(proband_alignment_classifications)
            # Log the results separately for every clade
            save_multiple_alignment_classifications(
                results_csv_file=results_csv_file,
                alignment_classifications=alignment_classifications,
                proband_number=initial_directory_name
            )
    return alignment_classifications


def edge_cut_maximum_alignable_subclade():
    perfect_data_dir_path = get_directory_path("Specify the path to the super directory with the simulated data:")
    perfect_data_dir_path = Path(perfect_data_dir_path)
    edge_cut_dir_path = get_directory_path("Specify the path to the directory with cut trees:")
    edge_cut_dir_path = Path(edge_cut_dir_path)
    os.chdir(edge_cut_dir_path)
    simulation_name = get_non_existing_path("Specify the simulation name:")
    simulation_path = edge_cut_dir_path / simulation_name
    os.makedirs(simulation_path, exist_ok=True)
    partial_results_filepath = simulation_path / partial_results_csv_filename
    # Iterate over the directories and identify the directories for alignment
    alignment_task_data = []
    for proband_number_directory in os.listdir(edge_cut_dir_path):
        # This is the directory that groups the data by the number of probands
        proband_number_directory_path = edge_cut_dir_path / proband_number_directory
        if not os.path.isdir(proband_number_directory_path):
            continue
        for original_tree_directory in os.listdir(proband_number_directory_path):
            # This is the directory storing all the data for the same original (simulated) tree
            original_tree_directory_path = proband_number_directory_path / original_tree_directory
            if not os.path.isdir(original_tree_directory_path):
                continue
            corresponding_perfect_data_path = (perfect_data_dir_path / proband_number_directory
                                               / original_tree_directory)
            if not os.path.isdir(corresponding_perfect_data_path):
                continue
            paths = get_paths_from_tree_pedigree_directory(corresponding_perfect_data_path)
            if not paths:
                continue
            pedigree_path, coalescent_tree_path = paths
            alignment_task_data.append((pedigree_path, coalescent_tree_path, original_tree_directory_path))
    # Run the alignments
    proband_number_to_results = defaultdict(list)
    with open(partial_results_filepath, "a") as partial_results_file:
        partial_results_file.write(csv_header)
        for alignment_task_data_entry in alignment_task_data:
            pedigree_path, coalescent_tree_path, original_tree_directory_path = alignment_task_data_entry
            # We are going to align the pedigree with multiple trees, so it's better to process the initial pedigree
            # once instead of processing a slightly smaller version for every tree
            # This, of course, cannot be efficient if the pedigree hasn't been reduced to the ascending
            # pedigree for the initial error-tree (which we use for edge cutting)
            pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path,
                                                                                 preprocess_graph=True)
            coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=coalescent_tree_path)
            root_vertex = coalescent_tree.get_root_vertex()
            del coalescent_tree
            root_vertex_information = RootVertexSpouses.get_root_vertex_information_from_pedigree(
                pedigree=pedigree, root_vertex_ploid_id=root_vertex
            )
            for modified_tree_directory in os.listdir(original_tree_directory_path):
                # This is the directory storing the data for the same modified tree.
                # Specifically, it stores all the cut trees derived from the modified tree
                modified_tree_directory_path = original_tree_directory_path / modified_tree_directory
                modified_tree_results_path = modified_tree_directory_path / simulation_name
                if not os.path.isdir(modified_tree_directory_path):
                    continue
                file_access = AlignmentResultsFileAccess(simulation_directory_path=modified_tree_results_path)
                task = MaximumAlignableSubcladeTask(
                    pedigree=pedigree,
                    tree_cut_clades_dir_path=modified_tree_directory_path,
                    root_vertex_info=root_vertex_information,
                    file_access=file_access,
                    simulation_subdirectory_path=modified_tree_results_path
                )
                proband_number, alignment_classification = task.run()
                if not proband_number:
                    continue
                csv_row = alignment_classification.csv_row(proband_number=str(proband_number))
                partial_results_file.write(csv_row)
                partial_results_file.flush()
                proband_number_to_results[proband_number].append(alignment_classification)
    results_csv_filepath = simulation_path / results_csv_filename
    global_alignment_classification = save_alignment_results_sorted_by_proband_number(
                                            results_csv_filepath=results_csv_filepath,
                                            proband_number_to_result=proband_number_to_results
    )
    os.remove(partial_results_filepath)
    return global_alignment_classification


def run_interactive_session():
    script_menu = ("Choose the running mode:\n"
                   "1) Specify the path to the initial tree-pedigree pair and a parent directory of error-simulated "
                   "pedigree directories\n"
                   "2) Specify the path to multiple tree-pedigree pair directories and a parent directory of "
                   "error-simulated pedigree directories\n"
                   "3) Specify the path to a parent directory containing tree-pedigree pair directories and "
                   "a parent directory of error-simulated pedigree directories\n"
                   "4) Specify the path to initial tree-pedigree pair and a parent directory of error-simulated "
                   "trees\n"
                   "5) Run the alignments for multiple tree-pedigree pairs and their corresponding error trees\n"
                   "6) (Ploid number simulation) Specify a super-parent directory whose subdirectories are "
                   "grouped by the proband number. Then, specify a super-parent directory with modified trees "
                   "directories (e.g. obtained from the create/resolve polytomy script)\n"
                   "7) (Edge cut simulation) Align the maximum alignable subclade\n"
                   "8) Run the alignments for a directory with multiple tree-pedigree subdirectories\n"
                   )
    menu_option = get_natural_number_input_in_bounds(script_menu, 1, 8)
    if menu_option < 7:
        prompt_top_founders()
    match menu_option:
        case 1:
            single_initial_data_directory_pedigree_error_directory_alignment_mode()
        case 2:
            specific_error_pedigree_directories_mode()
        case 3:
            children_tree_directories_error_pedigree_directory_mode()
        case 4:
            single_initial_data_directory_tree_error_directory_alignment_mode()
        case 5:
            multiple_tree_pedigree_directories_and_corresponding_tree_error_directories()
        case 6:
            create_or_resolve_polytomy_script()
        case 7:
            edge_cut_maximum_alignable_subclade()
        case 8:
            tree_pedigree_subdirectories()


def validate_required_arguments(arguments, required_args):
    missing_args = [arg for arg in required_args if getattr(arguments,
                                                            arg.replace('--', '').replace('-', '_')) is None]
    if missing_args:
        print(f"Missing required arguments for mode {arguments.mode}: {', '.join(missing_args)}")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description=script_help_description)
    parser.add_argument("--error-dir", type=str, required=True,
                        help="Path to the parent directory of error-simulated pedigree directories.")
    parser.add_argument("--mode", choices=["1", "2", "3"], required=True,
                        help="Specify the mode of operation. Choose from 1, 2, or 3.")
    parser.add_argument("--simulation-name", type=str, required=True,
                        help="Specify the simulation subpath.")

    # Mode 1 specific arguments
    parser.add_argument("--tree-path", type=str, help="Path to the coalescent tree (Mode 1)")
    parser.add_argument("--pedigree-path", type=str, help="Path to the coalescent tree (Mode 1)")

    # Mode 2 specific arguments
    parser.add_argument("--initial-dirs", nargs='+', help="Path to multiple tree-pedigree "
                                                          "pair directories (Mode 2)")

    # Mode 3 specific arguments
    parser.add_argument("--initial-parent-dir", type=str, help="Path to the parent directory containing "
                                                               "tree-pedigree pair directories (Mode 3)")
    parser.add_argument("--parallel", action='store_true', help="Enable parallel processing (Mode 3)")
    parser.add_argument("--max-workers", type=int, help="Maximum number of workers for "
                                                        "parallel processing (Mode 3)")

    # Parse the arguments
    arguments = parser.parse_args()
    # Define required arguments for each mode
    mode_args = {
        "1": ["--tree-path", "--pedigree-path"],
        "2": ["--initial-dirs"],
        "3": ["--initial-parent-dir", "--parallel", "--max-workers"]
    }
    shared_arguments = ["--error-dir", "--mode", "--simulation-name"]
    # Get the required arguments for the specified mode
    mode_specific_args = mode_args[arguments.mode]

    # Get all the arguments passed by the user
    passed_args = [arg for arg in vars(arguments) if getattr(arguments, arg) is not None]
    # Normalize argument names (convert underscores to hyphens for comparison)
    normalized_passed_args = [f"--{arg.replace('_', '-')}" for arg in passed_args]

    # Validate that only the required arguments for the selected mode are passed
    for arg in normalized_passed_args:
        if arg not in shared_arguments and arg not in mode_specific_args:
            print(f"Mode {arguments.mode} does not accept the argument: {arg}")
            sys.exit(1)

    # Check that the error directory exists and the simulation subpath is not taken (i.e. there is no folder with that
    # name under the error directory)
    error_directory_path = arguments.error_dir
    simulation_name = arguments.simulation_name
    if not utility.verify_folder_path(error_directory_path):
        print(f"The specified directory path {error_directory_path} either does not exist or is not a directory")
        return
    simulation_folder_path = Path(error_directory_path) / simulation_name
    if os.path.exists(simulation_folder_path):
        print(f"The specified simulation path {simulation_folder_path} already exists")
        return
    # Check the mode and print the corresponding information
    print(f"The specified mode is: {arguments.mode}")
    print(f"Error directory is: {error_directory_path}")
    selected_mode = arguments.mode
    match selected_mode:
        case "1":
            pedigree_path = arguments.pedigree_path
            tree_path = arguments.tree_path
            if tree_path and pedigree_path:
                for path in (tree_path, pedigree_path):
                    if not utility.verify_filepath(path):
                        print(f"Path {path} does not exist or is not a file")
                        return
                print(f"Running mode 1 with tree_path: {tree_path} and pedigree_path: {pedigree_path}")
                run_single_data_pedigree_error_directory_alignment(tree_path=tree_path,
                                                                   pedigree_no_errors_path=pedigree_path,
                                                                   error_pedigrees_folder_path=error_directory_path,
                                                                   input_simulation_subpath=simulation_name)
            else:
                print("Mode 1 requires both --tree-path and --pedigree-path.")
                return
        case "2":
            if arguments.initial_dirs:
                pedigree_tree_directory_paths = arguments.initial_dirs
                for directory in pedigree_tree_directory_paths:
                    if not utility.verify_folder_path(directory):
                        print(f"Path {directory} does not exist or is not a directory")
                        return
                print(f"Running mode 2 selected with initial_dirs: {arguments.initial_dirs}")
                run_specific_error_pedigree_directories(paths=pedigree_tree_directory_paths,
                                                        error_pedigrees_folder_path=error_directory_path,
                                                        simulation_subpath=simulation_name)
            else:
                print("Mode 2 requires --initial-dirs.")
                return
        case "3":
            if arguments.initial_parent_dir:
                parent_paths_directory = arguments.initial_parent_dir
                if not utility.verify_folder_path(parent_paths_directory):
                    print(f"Path {parent_paths_directory} does not exist or is not a directory")
                    return
                pedigree_tree_directory_paths = get_absolute_paths_to_subfolders(parent_paths_directory)
                if not pedigree_tree_directory_paths:
                    print(f"The directory {parent_paths_directory} does not contain any subdirectories")
                    return
                parallel_processing = arguments.parallel
                max_workers = arguments.max_workers
                if parallel_processing:
                    if max_workers is None or max_workers <= 1:
                        print("When --parallel is set, you must provide --max-workers with a value greater than 1.")
                        sys.exit(1)
                print(f"Running mode 3 selected with initial_parent_dir: {parent_paths_directory}")
                run_specific_error_pedigree_directories(paths=pedigree_tree_directory_paths,
                                                        error_pedigrees_folder_path=error_directory_path,
                                                        simulation_subpath=simulation_name,
                                                        parallelize=parallel_processing,
                                                        max_workers=max_workers)
            else:
                print("Mode 3 requires --initial-parent-dir.")
                return


def main():
    if len(sys.argv) == 1:  # Only the script name is present
        run_interactive_session()
    else:
        parse_arguments()


if __name__ == "__main__":
    main()
