from __future__ import annotations

import argparse
import os.path
import sys
import traceback
import warnings
from concurrent.futures import as_completed, ProcessPoolExecutor
from io import UnsupportedOperation

from alignment.graph_matcher import *
from scripts import utility
from scripts.run_alignment import save_alignment_result_to_files
from scripts.utility import *
from abc import ABC, abstractmethod

log_directory = "log_dir"
initial_alignment_dir_name = "initial_alignment"
general_result_filename = "simulation_result.txt"
detailed_result_filename = "detailed_simulation_result.txt"

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


class ErrorPedigreeAlignmentClassification:

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
        # The next fields keep track of the initial and new assignments to the root vertex
        self.total_number_of_root_vertex_individual_assignments = 0
        self.initial_root_vertex_individual_assignments_number = 0

    def add_results(self, other_result: ErrorPedigreeAlignmentClassification):
        self.simulations_with_solutions += other_result.simulations_with_solutions
        self.no_solutions_counter += other_result.no_solutions_counter
        self.only_individual_and_spouses_counter += other_result.only_individual_and_spouses_counter
        self.individual_and_non_spouse_present_counter += other_result.individual_and_non_spouse_present_counter
        self.individual_not_present_spouse_present_counter += (
            other_result.individual_not_present_spouse_present_counter)
        self.neither_individual_nor_spouse_present_counter += (
            other_result.neither_individual_nor_spouse_present_counter)
        self.total_number_of_root_vertex_individual_assignments += (
            other_result.total_number_of_root_vertex_individual_assignments)
        self.total_simulation_number += other_result.total_simulation_number


class ErrorAlignmentResultsFilesAccess:
    def __init__(self, general_result_filepath: Path, detailed_result_filepath: Path):
        self.general_result_filepath = general_result_filepath
        self.detailed_result_filepath = detailed_result_filepath
        # A simple counter keeping track of how many alignments have been saved to these files
        self.simulation_step = 0
        self.detailed_result_file = None
        self.general_result_file = None

    def open_files(self):
        self.general_result_file = open(self.general_result_filepath, 'a')
        self.detailed_result_file = open(self.detailed_result_filepath, 'a')

    def log_solutions(self, alignments: list[dict]):
        # Saving the information about the alignments into the detailed file
        self.detailed_result_file.write("#####################################\n")
        self.detailed_result_file.write(f"Simulation {self.simulation_step}\n")
        self.detailed_result_file.write(f"Number of incorrect alignments in this simulation: {len(alignments)}\n")
        for incorrect_alignment in alignments:
            self.detailed_result_file.write("------------------------------------\n")
            self.detailed_result_file.write(f"{incorrect_alignment}\n")
            self.detailed_result_file.write("------------------------------------\n")
        self.detailed_result_file.write("#####################################\n")
        self.detailed_result_file.flush()
        # Saving the information about the alignments into the general file
        self.general_result_file.write("#####################################\n")
        self.general_result_file.write(f"Simulation step: {self.simulation_step}\n")
        self.general_result_file.write(
            f"Number of alignments in this simulation {len(alignments)}\n")
        self.simulation_step += 1

    def close_files(self):
        if self.general_result_file:
            self.general_result_file.close()
            self.general_result_file = None
        if self.detailed_result_file:
            self.detailed_result_file.close()
            self.detailed_result_file = None


class ErrorAlignmentTask:
    def __init__(self, pedigree_path: Path, coalescent_tree_path: Path,
                 files_access: ErrorAlignmentResultsFilesAccess,
                 root_vertex_individual_id: int,
                 root_vertex_pedigree_individual_and_spouses: [int],
                 coalescent_vertex_id: int = None,
                 error_results: ErrorPedigreeAlignmentClassification = None,
                 simulation_subdirectory_path: Path = None):
        self.pedigree_path = pedigree_path
        self.coalescent_tree_path = coalescent_tree_path
        self.files_access = files_access
        if not error_results:
            error_results = ErrorPedigreeAlignmentClassification()
        self.error_results = error_results
        self.simulation_subdirectory_path = simulation_subdirectory_path
        self.root_vertex_individual_id = root_vertex_individual_id
        self.root_vertex_pedigree_individual_and_spouses = root_vertex_pedigree_individual_and_spouses
        self.coalescent_vertex_id = coalescent_vertex_id

    def run(self):
        print(f"Parsing {self.pedigree_path} [{self.coalescent_tree_path}]")
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=self.coalescent_tree_path)
        if default_initial_matching_mode == InitialMatchingMode.PLOID:
            ascending_genealogy_probands = coalescent_tree.probands
        else:
            proband_individuals = {proband // 2 for proband in coalescent_tree.probands}
            probands_all_ploids = [
                ploid
                for proband_individual in proband_individuals
                for ploid in [2 * proband_individual, 2 * proband_individual + 1]
            ]
            ascending_genealogy_probands = probands_all_ploids
        pedigree = (PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=self.pedigree_path,
                                                                              preprocess_graph=True,
                                                                              probands=ascending_genealogy_probands
                                                                              )
                    )
        if self.simulation_subdirectory_path:
            os.mkdir(self.simulation_subdirectory_path)
            log_directory_path = str(self.simulation_subdirectory_path / logs_default_directory_name)
            logger = MatcherLogger(logs_directory_path=log_directory_path)
        else:
            logger = None
        matcher = GraphMatcher(coalescent_tree=coalescent_tree,
                               processed_graph=pedigree,
                               logger=logger)
        result = matcher.find_mapping()
        if not self.coalescent_vertex_id:
            number_of_clades = len(result)
            if number_of_clades != 1:
                raise UnsupportedOperation("The current version of the script supports trees with only one clade")
            root_vertex_id = next(iter(result))
        else:
            root_vertex_id = self.coalescent_vertex_id
        # root_vertex_individual = root_vertex // 2
        # root_vertex_children_all_ploids = [
        #     child_ploid
        #     for child in pedigree.children_map[root_vertex]
        #     for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
        # ]
        # individual_and_spouses = {p // 2 for child_ploid in root_vertex_children_all_ploids
        #                           for p in pedigree.parents_map.get(child_ploid, [])}
        if self.simulation_subdirectory_path:
            alignments_dir_path = str(self.simulation_subdirectory_path)
            save_alignment_result_to_files(alignment_result=result,
                                           coalescent_tree=coalescent_tree,
                                           pedigree=pedigree,
                                           directory_path=alignments_dir_path)
        resulting_alignments = [d for lst in result.values() for d in lst]
        self.error_results.total_simulation_number += 1
        # Calculate the number of root assignments for this simulation
        simulation_root_vertex_individual_assignments = {x[root_vertex_id] // 2 for x in
                                                         resulting_alignments}
        self.files_access.open_files()
        self.files_access.log_solutions(alignments=resulting_alignments)
        if len(resulting_alignments) > 0:
            self.files_access.general_result_file.write(f"Root vertex individual assignments "
                                                        f"({len(simulation_root_vertex_individual_assignments)}):"
                                                        f"{simulation_root_vertex_individual_assignments}\n")
            # Can close the files earlier to let other processes work with them
            self.files_access.close_files()
            self.error_results.total_number_of_root_vertex_individual_assignments += (
                len(simulation_root_vertex_individual_assignments))
            # Classify the simulation result
            if self.root_vertex_individual_id in simulation_root_vertex_individual_assignments:
                if simulation_root_vertex_individual_assignments.issubset(
                        self.root_vertex_pedigree_individual_and_spouses):
                    self.error_results.only_individual_and_spouses_counter += 1
                else:
                    self.error_results.individual_and_non_spouse_present_counter += 1
            else:
                if simulation_root_vertex_individual_assignments.intersection(
                        self.root_vertex_pedigree_individual_and_spouses):
                    self.error_results.individual_not_present_spouse_present_counter += 1
                else:
                    self.error_results.neither_individual_nor_spouse_present_counter += 1
        else:
            self.error_results.no_solutions_counter += 1
        self.files_access.close_files()


def run_initial_alignment(coalescent_tree: CoalescentTree, initial_pedigree: PotentialMrcaProcessedGraph,
                          result_filepath: str, save_alignments_to_files: bool):
    result_filepath = Path(result_filepath)
    os.mkdir(result_filepath)
    if save_alignments_to_files:
        log_directory_path = result_filepath / logs_default_directory_name
        initial_logger = MatcherLogger(logs_directory_path=log_directory_path)
    else:
        initial_logger = None
    initial_matcher = GraphMatcher(coalescent_tree=coalescent_tree, processed_graph=initial_pedigree,
                                   logger=initial_logger)
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
        self.error_comparison_results = ErrorPedigreeAlignmentClassification()
        self.initial_coalescent_tree_path = initial_coalescent_tree_path
        self.initial_pedigree_path = initial_pedigree_path
        self.save_alignments_to_files = save_alignments_to_files
        self.simulation_folder_path = Path(simulation_folder_path)
        self.simulation_folder_path.mkdir(parents=True, exist_ok=True)
        detailed_result_filepath = self.simulation_folder_path / detailed_result_filename
        general_result_filepath = self.simulation_folder_path / general_result_filename
        self.file_access = ErrorAlignmentResultsFilesAccess(
            general_result_filepath=general_result_filepath,
            detailed_result_filepath=detailed_result_filepath
        )
        self.log_graph_paths()

    def log_graph_paths(self):
        general_result_file = open(self.file_access.general_result_filepath, 'a')
        general_result_file.write(f"The path to the initial, error-free pedigree is:"
                                  f" {self.initial_pedigree_path}\n")
        general_result_file.write(f"The path to the coalescent tree is: {self.initial_coalescent_tree_path}\n")
        general_result_file.write("#####################################\n")
        general_result_file.close()

    def log_overall_results(self):
        self.error_comparison_results.simulations_with_solutions = (
                self.error_comparison_results.total_simulation_number -
                self.error_comparison_results.no_solutions_counter)
        general_result_file = open(self.file_access.general_result_filepath, 'a')
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
        general_result_file.write(f"The number of simulations without the individual and all the spouses: "
                                  f"{self.error_comparison_results.neither_individual_nor_spouse_present_counter}\n")
        general_result_file.write("----------------------------------------------------------------------\n")
        general_result_file.write(f"Initial number of individual candidates for the root: "
                                  f"{self.error_comparison_results.initial_root_vertex_individual_assignments_number}\n")
        if self.error_comparison_results.simulations_with_solutions != 0:
            general_result_file.write(
                f"Average number of individual candidates for the root per simulation with solutions: "
                f"{self.error_comparison_results.total_number_of_root_vertex_individual_assignments / self.error_comparison_results.simulations_with_solutions}\n")
        general_result_file.close()

    def run_initial_alignment_and_log_general_graph_information(self) -> (int, [int]):
        initial_alignment_result_path = self.simulation_folder_path / initial_alignment_dir_name
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=self.initial_coalescent_tree_path)
        coalescent_tree.remove_unary_nodes()
        number_of_probands = len(coalescent_tree.probands)
        non_proband_vertices = len(coalescent_tree.vertex_to_level_map) - number_of_probands
        general_result_file = open(self.file_access.general_result_filepath, 'a')
        general_result_file.write(f"The number of probands is {number_of_probands}\n")
        general_result_file.write(f"The number of non-proband vertices is: {non_proband_vertices}\n")
        general_result_file.write("----------------------------------------------------------------------\n")
        if default_initial_matching_mode == InitialMatchingMode.PLOID:
            ascending_genealogy_probands = coalescent_tree.probands
        else:
            proband_individuals = {proband // 2 for proband in coalescent_tree.probands}
            probands_all_ploids = [
                ploid
                for proband_individual in proband_individuals
                for ploid in [2 * proband_individual, 2 * proband_individual + 1]
            ]
            ascending_genealogy_probands = probands_all_ploids
        initial_pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(
            filepath=self.initial_pedigree_path,
            probands=ascending_genealogy_probands,
            preprocess_graph=True
        )
        initial_result_dict = run_initial_alignment(coalescent_tree=coalescent_tree,
                                                    initial_pedigree=initial_pedigree,
                                                    save_alignments_to_files=self.save_alignments_to_files,
                                                    result_filepath=initial_alignment_result_path)
        number_of_clades = len(initial_result_dict)
        if number_of_clades != 1:
            raise UnsupportedOperation("The current version of the script supports trees with only one clade")
        root_vertex = next(iter(initial_result_dict))
        root_vertex_individual = root_vertex // 2
        root_vertex_children_all_ploids = [
            child_ploid
            for child in initial_pedigree.children_map[root_vertex]
            for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
        ]
        # In the phased case, we cannot guarantee that the other ploid of every root vertex
        # child ploid is present in the graph
        individual_and_spouses = {p // 2 for child_ploid in root_vertex_children_all_ploids
                                  for p in initial_pedigree.parents_map.get(child_ploid, [])}
        assert individual_and_spouses, "Invalid pedigree"
        valid_alignments = [d for lst in initial_result_dict.values() for d in lst]
        # Log the general information about the graphs and the initial alignment
        general_result_file.write(f"Total number of initial alignments is: {len(valid_alignments)}\n")
        # Calculating the number of individuals to which the clade's root is assigned
        initial_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in valid_alignments}
        self.error_comparison_results.initial_root_vertex_individual_assignments_number = (
            len(initial_root_vertex_individual_assignments)
        )
        general_result_file.write(f"Number of individual assignments for the root: "
                                  f"{self.error_comparison_results.initial_root_vertex_individual_assignments_number}\n"
                                  )
        general_result_file.write(f"Clade root individual: {root_vertex_individual}\n")
        general_result_file.write(f"Clade root and their spouses: {individual_and_spouses}\n")
        general_result_file.write(
            f"Individual assignments for the root: {initial_root_vertex_individual_assignments}\n")
        general_result_file.flush()
        general_result_file.close()
        return root_vertex, root_vertex_individual, individual_and_spouses

    def run_alignments(self) -> ErrorPedigreeAlignmentClassification:
        alignment_tasks: [ErrorAlignmentTask] = self.get_alignment_tasks()
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

    def get_alignment_tasks(self) -> [ErrorAlignmentTask]:
        root_vertex_id, root_vertex_individual_id, individual_and_spouses = (
            self.run_initial_alignment_and_log_general_graph_information())
        current_directory = Path(self.error_pedigrees_folder)
        alignment_tasks = []
        for subdirectory in os.listdir(self.error_pedigrees_folder):
            simulation_subdirectory_path = self.simulation_folder_path / subdirectory
            subdirectory_path = current_directory / subdirectory
            pedigree_file = get_unique_filename_with_specified_extension(directory_path=subdirectory_path,
                                                                         extension=".pedigree")

            pedigree_filepath = subdirectory_path / pedigree_file
            alignment_task = ErrorAlignmentTask(
                pedigree_path=pedigree_filepath,
                coalescent_tree_path=self.initial_coalescent_tree_path,
                files_access=self.file_access,
                simulation_subdirectory_path=simulation_subdirectory_path,
                error_results=self.error_comparison_results,
                root_vertex_individual_id=root_vertex_individual_id,
                root_vertex_pedigree_individual_and_spouses=individual_and_spouses
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

    def get_alignment_tasks(self) -> [ErrorAlignmentTask]:
        root_vertex_id, root_vertex_individual_id, individual_and_spouses = (
            self.run_initial_alignment_and_log_general_graph_information())
        # Loop over the error trees and generate the alignment tasks
        current_directory = Path(self.tree_error_directory_path)
        alignment_tasks = []
        for coalescent_tree_file in os.listdir(self.tree_error_directory_path):
            coalescent_tree_filepath = current_directory / coalescent_tree_file
            if not os.path.isfile(coalescent_tree_filepath):
                continue
            simulation_subdirectory_path = self.simulation_folder_path / coalescent_tree_file
            alignment_task = ErrorAlignmentTask(
                pedigree_path=self.initial_pedigree_path,
                coalescent_tree_path=coalescent_tree_filepath,
                files_access=self.file_access,
                simulation_subdirectory_path=simulation_subdirectory_path,
                error_results=self.error_comparison_results,
                root_vertex_individual_id=root_vertex_individual_id,
                root_vertex_pedigree_individual_and_spouses=individual_and_spouses,
                coalescent_vertex_id=root_vertex_id
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
    error_pedigree_alignment_comparison.run_alignments()


def run_single_data_tree_error_directory_alignment(tree_path: str, pedigree_path: str,
                                                   error_trees_folder_path: str, simulation_subpath: str) \
        -> ErrorPedigreeAlignmentClassification:
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
                 error_pedigrees_folder_path: str | Path, simulation_subpath: str) -> ErrorPedigreesAlignmentComparison:
    absolute_path_to_process = simulation_folder_path / path_to_process
    paths = get_paths_from_tree_pedigree_directory(absolute_path_to_process)
    if not paths:
        warnings.warn(f"The {path_to_process} directory is invalid. Skipping.")
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


def process_path_and_run_initial_alignments(path_to_process: str | Path, simulation_folder_path: str | Path,
                                            error_pedigrees_folder_path: str | Path, simulation_name: str
                                            ) -> [ErrorAlignmentTask]:
    tree_alignments: ErrorPedigreesAlignmentComparison = process_path(path_to_process, simulation_folder_path,
                                                                      error_pedigrees_folder_path, simulation_name)
    alignment_tasks: [ErrorAlignmentTask] = tree_alignments.get_alignment_tasks()
    return tree_alignments, alignment_tasks


def save_error_alignments_to_file(filepath: str, error_alignments_results: ErrorPedigreeAlignmentClassification):
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
        results_file.write(f"The number of simulations without the individual and all the spouses: "
                           f"{error_alignments_results.neither_individual_nor_spouse_present_counter}\n")
        results_file.close()


def run_specific_error_pedigree_directories_iterative(paths: list[str], error_pedigrees_folder_path: str,
                                                      simulation_subpath: str):
    error_pedigree_parent_directory = Path(error_pedigrees_folder_path).parent
    simulation_folder_path = error_pedigree_parent_directory / simulation_subpath
    os.makedirs(simulation_folder_path, exist_ok=True)
    total_results = ErrorPedigreeAlignmentClassification()
    for path_to_process in paths:
        error_alignment = process_path(path_to_process=path_to_process,
                                       simulation_folder_path=simulation_folder_path,
                                       error_pedigrees_folder_path=error_pedigrees_folder_path,
                                       simulation_subpath=simulation_subpath)
        tree_results = error_alignment.run_alignments()
        total_results.add_results(tree_results)

    total_results_filepath = simulation_folder_path / "total_results.txt"
    save_error_alignments_to_file(total_results_filepath, total_results)
    return total_results


def get_alignments_tasks(error_alignment_comparison: ErrorPedigreesAlignmentComparison):
    return error_alignment_comparison, error_alignment_comparison.get_alignment_tasks()


def process_alignment_task(path_alignment_task: ErrorAlignmentTask):
    path_alignment_task.run()
    return path_alignment_task


def run_specific_error_pedigree_directories_parallel(paths: list[str], error_pedigrees_folder_path: str,
                                                     simulation_subpath: str, max_workers: int):
    error_pedigree_parent_directory = Path(error_pedigrees_folder_path).parent
    simulation_directory_path = error_pedigree_parent_directory / simulation_subpath
    os.makedirs(simulation_directory_path)

    total_results = ErrorPedigreeAlignmentClassification()
    tree_alignment_comparisons: [ErrorPedigreesAlignmentComparison] = []
    # Separating all the pedigree-tree directories into separate comparison objects
    for path in paths:
        tree_alignment_comparison: ErrorPedigreesAlignmentComparison = process_path(path, simulation_directory_path,
                                                                                    error_pedigrees_folder_path,
                                                                                    simulation_subpath)
        tree_alignment_comparisons.append(tree_alignment_comparison)
    alignment_tasks: [ErrorAlignmentTask] = []
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
    for tree_alignment_comparison in tree_alignment_comparisons:
        tree_alignment_comparison.log_overall_results()
        total_results.add_results(tree_alignment_comparison.error_comparison_results)
    total_results_filepath = simulation_directory_path / "total_results.txt"
    save_error_alignments_to_file(total_results_filepath, total_results)
    return total_results


def run_specific_error_pedigree_directories(paths: list[str], error_pedigrees_folder_path: str,
                                            simulation_subpath: str, parallelize: bool = True,
                                            max_workers: int = 0):
    if parallelize:
        return run_specific_error_pedigree_directories_parallel(paths=paths,
                                                                error_pedigrees_folder_path=error_pedigrees_folder_path,
                                                                simulation_subpath=simulation_subpath,
                                                                max_workers=max_workers)
    return run_specific_error_pedigree_directories_iterative(paths=paths,
                                                             error_pedigrees_folder_path=error_pedigrees_folder_path,
                                                             simulation_subpath=simulation_subpath)


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
    tree_pedigree_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory to"
                                                             " tree-pedigree directories:")
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


def multiple_tree_pedigree_directories_and_corresponding_tree_error_directories():
    tree_pedigree_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory "
                                                             "containing tree-pedigree directories:")
    tree_pedigree_parent_directory_path = Path(tree_pedigree_parent_directory_path)
    error_trees_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory containing"
                                                           " error trees subdirectories:")
    error_trees_parent_directory_path = Path(error_trees_parent_directory_path)
    input_simulation_subpath = get_non_existing_path("Specify the simulation subpath:")
    input_simulation_subpath = Path(input_simulation_subpath)
    simulation_path = error_trees_parent_directory_path / input_simulation_subpath
    total_error_comparison_results = ErrorPedigreeAlignmentClassification()
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
        total_error_comparison_results.add_results(tree_pedigree_error_comparison)
    total_results_filepath = simulation_path / "total_results.txt"
    save_error_alignments_to_file(total_results_filepath, total_error_comparison_results)
    return total_error_comparison_results


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
                   )
    menu_option = get_natural_number_input_in_bounds(script_menu, 1, 5)
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
