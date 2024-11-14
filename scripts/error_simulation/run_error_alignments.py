from __future__ import annotations

import argparse
import os.path
import sys
import traceback
import warnings
from concurrent.futures import as_completed, ProcessPoolExecutor
from io import UnsupportedOperation
from threading import Lock

from alignment.graph_matcher import *
from scripts import utility
from scripts.run_alignment import save_alignment_result_to_files, get_alignments_proband_distance_probands_ignored
from scripts.utility import *

log_directory = "log_dir"
initial_alignment_dir_name = "initial_alignment"
general_result_filename = "simulation_result.txt"
detailed_result_filename = "detailed_simulation_result.txt"


class ErrorPedigreeAlignmentComparison:
    class ErrorPedigreeAlignmentClassification:

        def __init__(self):
            self.total_number_of_solutions = 0
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

        def add_results(self, other_result: ErrorPedigreeAlignmentComparison.ErrorPedigreeAlignmentClassification):
            self.total_number_of_solutions += other_result.total_number_of_solutions
            self.no_solutions_counter += other_result.no_solutions_counter
            self.only_individual_and_spouses_counter += other_result.only_individual_and_spouses_counter
            self.individual_and_non_spouse_present_counter += other_result.individual_and_non_spouse_present_counter
            self.individual_not_present_spouse_present_counter += other_result.individual_not_present_spouse_present_counter

    def __init__(self, coalescent_tree_path: str, error_free_pedigree_path: str,
                 error_pedigrees_folder: str, simulation_folder_subpath: str,
                 save_alignments_to_files: bool = True):
        self.coalescent_tree_path = coalescent_tree_path
        self.error_free_pedigree_path = error_free_pedigree_path
        self.error_pedigrees_folder = error_pedigrees_folder
        self.save_alignments_to_files = save_alignments_to_files
        self.coalescent_tree = None
        self.initial_pedigree = None
        self.initial_average_distance_to_identity = 0
        self.initial_root_vertex_individual_assignments_number = 0
        self.total_distance_to_solution = 0
        self.total_number_of_root_vertex_individual_assignments = 0
        self.new_individual_assignments_number = 0
        self.new_individual_simulations_number = 0
        self.error_pedigree_results = None
        error_pedigree_parent_directory = Path(self.error_pedigrees_folder).parent
        simulation_directory = error_pedigree_parent_directory / simulation_folder_subpath
        simulation_directory.mkdir(parents=True, exist_ok=True)
        self.simulation_folder_path = simulation_directory
        self.detailed_result_file = open(self.simulation_folder_path / detailed_result_filename, 'w')
        self.general_result_file = open(self.simulation_folder_path / general_result_filename, 'w')
        self.log_graph_paths()

    def log_graph_paths(self):
        self.general_result_file.write(f"The path to the initial, error-free pedigree is:"
                                       f" {self.error_free_pedigree_path}\n")
        self.general_result_file.write(f"The path to the coalescent tree is: {self.coalescent_tree_path}\n")
        self.general_result_file.write("#####################################\n")

    def log_solutions(self, alignments: list[dict], simulation_step: int):
        self.detailed_result_file.write("#####################################\n")
        self.detailed_result_file.write(f"Simulation {simulation_step}\n")
        self.detailed_result_file.write(f"Number of incorrect alignments in this simulation: {len(alignments)}\n")
        for incorrect_alignment in alignments:
            self.detailed_result_file.write("------------------------------------\n")
            self.detailed_result_file.write(f"{incorrect_alignment}\n")
            self.detailed_result_file.write("------------------------------------\n")
        self.detailed_result_file.write("#####################################\n")
        self.detailed_result_file.flush()

    def log_overall_results(self, number_of_simulation_steps: int):
        number_of_probands = len(self.coalescent_tree.probands)
        non_proband_vertices = len(self.coalescent_tree.vertex_to_level_map) - number_of_probands
        if self.error_pedigree_results.total_number_of_solutions != 0:
            alignments_to_identity_average_distance = (self.total_distance_to_solution
                                                       / self.error_pedigree_results.total_number_of_solutions)
            delta = alignments_to_identity_average_distance - self.initial_average_distance_to_identity

            non_proband_vertices = len(self.coalescent_tree.vertex_to_level_map) - number_of_probands
            self.general_result_file.write(
                f"Initial average distance to the identity is: "
                f"{self.initial_average_distance_to_identity}\n")
            self.general_result_file.write(f"Average distance from simulation solutions to the identity "
                                           f"is: {alignments_to_identity_average_distance}\n")
            self.general_result_file.write(f"Average distance delta is: {delta}\n")
            self.general_result_file.write(
                f"Average distance is {alignments_to_identity_average_distance / non_proband_vertices} "
                f"from the number of non-proband vertices in the coalescent tree\n")
        self.general_result_file.write(f"The number of probands is {number_of_probands}\n")
        self.general_result_file.write(f"The number of non-proband vertices is: {non_proband_vertices}\n")
        self.general_result_file.write("----------------------------------------------------------------------\n")
        self.error_pedigree_results.simulations_with_solutions = (number_of_simulation_steps -
                                                                  self.error_pedigree_results.no_solutions_counter)
        self.general_result_file.write(f"The number of simulations with solutions: "
                                       f"{self.error_pedigree_results.simulations_with_solutions}\n")
        self.general_result_file.write(f"The number of simulations with no solutions: "
                                       f"{self.error_pedigree_results.no_solutions_counter}\n")
        self.general_result_file.write(f"The number of simulations with only the individual and spouses: "
                                       f"{self.error_pedigree_results.only_individual_and_spouses_counter}\n")
        self.general_result_file.write(f"The number of simulations with the individual and a non-spouse present: "
                                       f"{self.error_pedigree_results.individual_and_non_spouse_present_counter}\n")
        self.general_result_file.write(f"The number of simulations without the individual, but with a spouse: "
                                       f"{self.error_pedigree_results.individual_not_present_spouse_present_counter}\n")
        self.general_result_file.write(f"The number of simulations without the individual and all the spouses: "
                                       f"{self.error_pedigree_results.neither_individual_nor_spouse_present_counter}\n")
        self.general_result_file.write("----------------------------------------------------------------------\n")
        self.general_result_file.write(f"Initial number of individual candidates for the root: "
                                       f"{self.initial_root_vertex_individual_assignments_number}\n")
        if self.error_pedigree_results.simulations_with_solutions != 0:
            self.general_result_file.write(
                f"Average number of individual candidates for the root per simulation with solutions: "
                f"{self.total_number_of_root_vertex_individual_assignments / self.error_pedigree_results.simulations_with_solutions}\n")
            self.general_result_file.write(f"New individual root assignments: "
                                           f"{self.new_individual_assignments_number}\n")
            self.general_result_file.write(
                f"Average number of new individual root assignments per simulation with solutions: "
                f"{self.new_individual_assignments_number / self.error_pedigree_results.simulations_with_solutions}\n")
            self.general_result_file.write(f"Number of new individual root assignments simulations: "
                                           f"{self.new_individual_simulations_number}\n")

    def run_initial_alignment(self):
        initial_alignment_path = self.simulation_folder_path / initial_alignment_dir_name
        os.mkdir(initial_alignment_path)
        if self.save_alignments_to_files:
            log_directory_path = initial_alignment_path / logs_default_directory_name
            initial_logger = MatcherLogger(logs_directory_path=log_directory_path)
        else:
            initial_logger = None
        initial_matcher = GraphMatcher(coalescent_tree=self.coalescent_tree, processed_graph=self.initial_pedigree,
                                       logger=initial_logger)
        initial_result_dict = initial_matcher.find_mapping()
        if self.save_alignments_to_files:
            save_alignment_result_to_files(alignment_result=initial_result_dict,
                                           coalescent_tree=self.coalescent_tree,
                                           pedigree=self.initial_pedigree,
                                           directory_path=initial_alignment_path)
        return initial_result_dict

    def run_alignments(self) -> ErrorPedigreeAlignmentClassification:
        self.error_pedigree_results = ErrorPedigreeAlignmentComparison.ErrorPedigreeAlignmentClassification()
        self.coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=self.coalescent_tree_path)
        self.coalescent_tree.remove_unary_nodes()
        self.initial_pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(
            filepath=self.error_free_pedigree_path,
            initialize_levels=False,
            initialize_ancestor_maps=False
        )
        if default_initial_matching_mode == InitialMatchingMode.PLOID:
            ascending_genealogy_probands = self.coalescent_tree.probands
        else:
            proband_individuals = {proband // 2 for proband in self.coalescent_tree.probands}
            probands_all_ploids = [
                ploid
                for proband_individual in proband_individuals
                for ploid in [2 * proband_individual, 2 * proband_individual + 1]
            ]
            ascending_genealogy_probands = probands_all_ploids
        self.initial_pedigree.reduce_to_ascending_genealogy(probands=ascending_genealogy_probands,
                                                            recalculate_levels=True)
        self.initial_pedigree.initialize_potential_mrca_map()
        initial_result_dict = self.run_initial_alignment()
        number_of_clades = len(initial_result_dict)
        if number_of_clades != 1:
            raise UnsupportedOperation("The current version of the script supports trees with only one clade")
        root_vertex = next(iter(initial_result_dict))
        root_vertex_individual = root_vertex // 2
        root_vertex_children_all_ploids = [
            child_ploid
            for child in self.initial_pedigree.children_map[root_vertex]
            for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
        ]
        # In the phased case, we cannot guarantee that the other ploid of every root vertex
        # child ploid is present in the graph
        individual_and_spouses = {p // 2 for child_ploid in root_vertex_children_all_ploids
                                  for p in self.initial_pedigree.parents_map.get(child_ploid, [])}
        assert individual_and_spouses, "Zero assignments found for the root!"
        valid_alignments = [d for lst in initial_result_dict.values() for d in lst]
        # Letting the garbage collector to clean this object
        del initial_result_dict
        del self.initial_pedigree
        # Gathering the statistics about the original alignments
        identity_solution = self.coalescent_tree.get_identity_solution()
        assert identity_solution in valid_alignments
        initial_average_distance_to_identity = sum(
            get_alignments_proband_distance_probands_ignored(alignment, identity_solution,
                                                             self.coalescent_tree.probands)
            for alignment in valid_alignments
        ) / len(valid_alignments)
        self.general_result_file.write(f"Total number of initial alignments is: {len(valid_alignments)}\n")
        self.general_result_file.write(f"Average distance to the identity is: {initial_average_distance_to_identity}\n")
        # Calculating the number of individuals to which the clade's root is assigned
        initial_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in valid_alignments}
        self.initial_root_vertex_individual_assignments_number += (
            len(initial_root_vertex_individual_assignments))
        self.general_result_file.write(f"Number of individual assignments for the root:"
                                       f" {self.initial_root_vertex_individual_assignments_number}\n")
        self.general_result_file.write(f"Clade root individual: {root_vertex_individual}\n")
        self.general_result_file.write(f"Clade root and their spouses: {individual_and_spouses}\n")
        self.general_result_file.write(
            f"Individual assignments for the root: {initial_root_vertex_individual_assignments}\n")
        self.general_result_file.flush()

        # Run the alignments on the pedigrees with errors
        current_directory = Path(self.error_pedigrees_folder)
        print(f"Starting running the alignments for {self.coalescent_tree_path}")
        simulation_counter = 0
        try:
            for subdirectory in os.listdir(self.error_pedigrees_folder):
                subdirectory_path = current_directory / subdirectory
                simulation_counter += 1
                for file in os.listdir(subdirectory_path):
                    if file.endswith(".pedigree"):
                        filepath = subdirectory_path / file
                        print(f"Parsing {filepath} [{self.coalescent_tree_path}]")
                        potential_mrca_graph = (PotentialMrcaProcessedGraph.
                                                get_processed_graph_from_file(filepath=filepath,
                                                                              initialize_levels=False,
                                                                              initialize_ancestor_maps=False))
                        print("Reducing to the ascending genealogy")
                        potential_mrca_graph.reduce_to_ascending_genealogy(probands=ascending_genealogy_probands,
                                                                           recalculate_levels=True)
                        potential_mrca_graph.initialize_potential_mrca_map()
                        simulation_subdirectory_path = self.simulation_folder_path / subdirectory
                        if self.save_alignments_to_files:
                            os.mkdir(simulation_subdirectory_path)
                            log_directory_path = str(simulation_subdirectory_path / logs_default_directory_name)
                            logger = MatcherLogger(logs_directory_path=log_directory_path)
                        else:
                            logger = None
                        matcher = GraphMatcher(coalescent_tree=self.coalescent_tree,
                                               processed_graph=potential_mrca_graph,
                                               logger=logger)
                        result = matcher.find_mapping()
                        if self.save_alignments_to_files:
                            alignments_dir_path = str(simulation_subdirectory_path)
                            save_alignment_result_to_files(alignment_result=result,
                                                           coalescent_tree=self.coalescent_tree,
                                                           pedigree=potential_mrca_graph,
                                                           directory_path=alignments_dir_path)
                        resulting_alignments = [d for lst in result.values() for d in lst]
                        self.log_solutions(resulting_alignments, simulation_counter)
                        self.general_result_file.write("#####################################\n")
                        self.general_result_file.write(f"Simulation step: {simulation_counter}\n")
                        self.general_result_file.write(
                            f"Number of alignments in this simulation {len(resulting_alignments)}\n")
                        if len(resulting_alignments) > 0:
                            self.error_pedigree_results.total_number_of_solutions += len(resulting_alignments)
                            self.total_distance_to_solution += sum(
                                get_alignments_proband_distance_probands_ignored(alignment, identity_solution,
                                                                                 self.coalescent_tree.probands)
                                for alignment in resulting_alignments
                            )
                            (self.general_result_file.write(
                                f"New average distance from simulation solutions to the"
                                f" identity is: "
                                f"{self.total_distance_to_solution / self.error_pedigree_results.total_number_of_solutions}\n")
                            )
                            # Calculate the number of root assignments for this simulation
                            simulation_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in
                                                                             resulting_alignments}
                            self.general_result_file.write(f"Root vertex individual assignments "
                                                           f"({len(simulation_root_vertex_individual_assignments)}): "
                                                           f"{simulation_root_vertex_individual_assignments}\n")
                            self.total_number_of_root_vertex_individual_assignments += (
                                len(simulation_root_vertex_individual_assignments))
                            self.new_individual_assignments_number += len(
                                {x for x in simulation_root_vertex_individual_assignments
                                 if x not in initial_root_vertex_individual_assignments})
                            if self.new_individual_assignments_number > 0:
                                self.new_individual_simulations_number += 1
                            # Classify the simulation result
                            if root_vertex_individual in simulation_root_vertex_individual_assignments:
                                if simulation_root_vertex_individual_assignments.issubset(individual_and_spouses):
                                    self.error_pedigree_results.only_individual_and_spouses_counter += 1
                                else:
                                    self.error_pedigree_results.individual_and_non_spouse_present_counter += 1
                            else:
                                if simulation_root_vertex_individual_assignments.intersection(individual_and_spouses):
                                    self.error_pedigree_results.individual_not_present_spouse_present_counter += 1
                                else:
                                    self.error_pedigree_results.neither_individual_nor_spouse_present_counter += 1
                        else:
                            self.error_pedigree_results.no_solutions_counter += 1
                        self.general_result_file.write("#####################################\n")
                        self.general_result_file.flush()
                        break
        except KeyboardInterrupt:
            print("Stop the simulation, log the final results")
        except Exception as ex:
            print(f"Exception occurred: {ex}")
            traceback.print_exc()

        self.log_overall_results(number_of_simulation_steps=simulation_counter)
        self.detailed_result_file.close()
        self.general_result_file.close()
        return self.error_pedigree_results


def run_single_parent_error_directory_alignment(tree_path: str, pedigree_no_errors_path: str,
                                                error_pedigrees_folder_path: str, input_simulation_name: str):
    os.chdir(error_pedigrees_folder_path)
    os.chdir("..")
    error_pedigree_alignment_comparison = ErrorPedigreeAlignmentComparison(
        coalescent_tree_path=tree_path,
        error_pedigrees_folder=error_pedigrees_folder_path,
        error_free_pedigree_path=pedigree_no_errors_path,
        simulation_folder_subpath=input_simulation_name)
    error_pedigree_alignment_comparison.run_alignments()


def single_tree_parent_error_directory_alignment_mode():
    current_working_directory = Path.cwd()
    tree_path = get_file_path("Specify the absolute path to the coalescent tree:")
    pedigree_no_errors_path = get_file_path("Specify the absolute path to the initial, error-free pedigree:")
    error_pedigrees_folder_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    os.chdir(error_pedigrees_folder_path)
    input_simulation_name = get_non_existing_directory_path("Specify the simulation name:")
    os.chdir(current_working_directory)
    run_single_parent_error_directory_alignment(tree_path=tree_path, pedigree_no_errors_path=pedigree_no_errors_path,
                                                error_pedigrees_folder_path=error_pedigrees_folder_path,
                                                input_simulation_name=input_simulation_name)


def get_unique_folder(base_name: str, directory: Path = Path.cwd()) -> Path:
    # Initial attempt with the base name
    unique_name = base_name
    count = 1

    # Check if the directory exists, and add a prefix if it does
    while (directory / unique_name).exists():
        unique_name = f"{count}_{base_name}"
        count += 1
    return unique_name


def process_path(path_to_process, simulation_folder_path, error_pedigrees_folder_path, simulation_name):
    absolute_path_to_process = simulation_folder_path / path_to_process
    files = list(absolute_path_to_process.glob("*"))
    files = [file.resolve() for file in files if file.is_file()]

    if len(files) == 2:
        pedigree_path = next((file for file in files if file.suffix == ".pedigree"), None)
        tree_path = next((file for file in files if file != pedigree_path), None)

        if not pedigree_path or not tree_path:
            warnings.warn("Either the pedigree or the tree file is missing, skipping the simulation")
            return None
    else:
        warnings.warn(f"The {path_to_process} directory does not contain exactly two files. Skipping.")
        return None

    simulation_folder_basename = os.path.basename(path_to_process)
    simulation_folder_unique_name = get_unique_folder(base_name=simulation_folder_basename,
                                                      directory=simulation_folder_path)
    simulation_path = str(Path(simulation_name) / Path(simulation_folder_unique_name))
    error_pedigree_alignment_comparison = ErrorPedigreeAlignmentComparison(
        coalescent_tree_path=tree_path,
        error_pedigrees_folder=error_pedigrees_folder_path,
        error_free_pedigree_path=pedigree_path,
        simulation_folder_subpath=simulation_path
    )
    return error_pedigree_alignment_comparison.run_alignments()


def save_error_alignments_to_file(filepath: str,
                                  error_alignments_results: ErrorPedigreeAlignmentComparison.ErrorPedigreeAlignmentClassification):
    with open(filepath, "w") as results_file:
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
                                                      simulation_name: str):
    error_pedigree_parent_directory = Path(error_pedigrees_folder_path).parent
    simulation_folder_path = error_pedigree_parent_directory / simulation_name
    os.mkdir(simulation_folder_path)
    total_results = ErrorPedigreeAlignmentComparison.ErrorPedigreeAlignmentClassification()
    for path_to_process in paths:
        tree_results = process_path(path_to_process=path_to_process,
                                    simulation_folder_path=simulation_folder_path,
                                    error_pedigrees_folder_path=error_pedigrees_folder_path,
                                    simulation_name=simulation_name)
        if tree_results:
            total_results.add_results(tree_results)

    total_results_filepath = simulation_folder_path / f"{simulation_name}.txt"
    save_error_alignments_to_file(total_results_filepath, total_results)
    return total_results


def run_specific_error_pedigree_directories_parallel(paths: list[str], error_pedigrees_folder_path: str,
                                                     simulation_name: str):
    error_pedigree_parent_directory = Path(error_pedigrees_folder_path).parent
    simulation_folder_path = error_pedigree_parent_directory / simulation_name
    os.mkdir(simulation_folder_path)

    total_results = ErrorPedigreeAlignmentComparison.ErrorPedigreeAlignmentClassification()
    lock = Lock()

    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(process_path, path, simulation_folder_path, error_pedigrees_folder_path, simulation_name)
            for path in paths]
        for future in as_completed(futures):
            tree_results = future.result()
            if tree_results:
                with lock:
                    total_results.add_results(tree_results)
    total_results_filepath = simulation_folder_path / f"{simulation_name}.txt"
    save_error_alignments_to_file(total_results_filepath, total_results)
    return total_results


def run_specific_error_pedigree_directories(paths: list[str], error_pedigrees_folder_path: str,
                                            simulation_name: str, parallelize: bool = True):
    if parallelize:
        return run_specific_error_pedigree_directories_parallel(paths, error_pedigrees_folder_path, simulation_name)
    return run_specific_error_pedigree_directories_iterative(paths, error_pedigrees_folder_path,
                                                             simulation_name)


def specific_error_pedigree_directories_mode():
    current_working_directory = Path.cwd()
    paths = get_folder_paths("Specify the absolute paths to all the tree-pedigree directories")
    error_pedigrees_folder_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    os.chdir(error_pedigrees_folder_path)
    input_simulation_name = get_non_existing_directory_path("Specify the simulation name:")
    os.chdir(current_working_directory)
    run_specific_error_pedigree_directories(paths=paths,
                                            error_pedigrees_folder_path=error_pedigrees_folder_path,
                                            simulation_name=input_simulation_name)


def get_absolute_paths_to_subfolders(directory_path: str) -> list[str]:
    return [
        os.path.abspath(os.path.join(directory_path, subdir))
        for subdir in os.listdir(directory_path)
        if os.path.isdir(os.path.join(directory_path, subdir))
    ]


def children_tree_directories_error_pedigree_directory_mode():
    current_working_directory = Path.cwd()
    error_pedigrees_folder_path = get_directory_path("Specify the absolute path to the error pedigrees directory:")
    tree_pedigree_parent_directory_path = get_directory_path("Specify the absolute path to a parent directory to"
                                                             " tree-pedigree directories:")
    os.chdir(error_pedigrees_folder_path)
    input_simulation_name = get_non_existing_directory_path("Specify the simulation name:")
    tree_pedigree_paths = get_absolute_paths_to_subfolders(tree_pedigree_parent_directory_path)
    os.chdir(current_working_directory)
    run_specific_error_pedigree_directories(paths=tree_pedigree_paths,
                                            error_pedigrees_folder_path=error_pedigrees_folder_path,
                                            simulation_name=input_simulation_name)


def run_interactive_session():
    script_menu = ("Choose the running mode:\n"
                   "1) Specify the path to the initial tree-pedigree pair and a parent directory of error-simulated "
                   "pedigree directories\n"
                   "2) Specify the path to multiple tree-pedigree pair directories and a parent directory of "
                   "error-simulated pedigree directories\n"
                   "3) Specify the path to a parent directory containing tree-pedigree pair directories and "
                   "a parent directory of error-simulated pedigree directories\n"
                   )
    menu_option = get_natural_number_input_in_bounds(script_menu, 1, 3)
    match menu_option:
        case 1:
            single_tree_parent_error_directory_alignment_mode()
        case 2:
            specific_error_pedigree_directories_mode()
        case 3:
            children_tree_directories_error_pedigree_directory_mode()


def parse_arguments():
    parser = argparse.ArgumentParser(description="Choose the running mode.")
    parser.add_argument("--error-dir", type=str, required=True,
                        help="Path to the parent directory of error-simulated pedigree directories.")
    parser.add_argument("--mode", choices=["1", "2", "3"], required=True,
                        help="Specify the mode of operation. Choose from 1, 2, or 3.")
    parser.add_argument("--simulation-name", type=str, required=True,
                        help="Specify the simulation name.")

    # Mode 1 specific arguments
    parser.add_argument("--tree-path", type=str, help="Path to the coalescent tree")
    parser.add_argument("--pedigree-path", type=str, help="Path to the coalescent tree")

    # Mode 2 specific arguments
    parser.add_argument("--initial-dirs", nargs='+', help="Path to multiple tree-pedigree pair directories")

    # Mode 3 specific arguments
    parser.add_argument("--initial-parent-dir", type=str, help="Path to the parent directory containing "
                                                               "tree-pedigree pair directories")

    # Parse the arguments
    arguments = parser.parse_args()
    # Define required arguments for each mode
    mode_args = {
        "1": ["--tree-path", "--pedigree-path"],
        "2": ["--initial-dirs"],
        "3": ["--initial-parent-dir"]
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

    # Check that the error directory exists and the simulation name is not taken (i.e. there is no folder with that
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
    if arguments.mode == "1":
        pedigree_path = arguments.pedigree_path
        tree_path = arguments.tree_path
        if tree_path and pedigree_path:
            for path in (tree_path, pedigree_path):
                if not utility.verify_filepath(path):
                    print(f"Path {path} does not exist or is not a file")
                    return
            print(f"Running mode 1 with tree_path: {tree_path} and pedigree_path: {pedigree_path}")
            run_single_parent_error_directory_alignment(tree_path=tree_path,
                                                        pedigree_no_errors_path=pedigree_path,
                                                        error_pedigrees_folder_path=error_directory_path,
                                                        input_simulation_name=simulation_name)
        else:
            print("Mode 1 requires both --tree-path and --pedigree-path.")
            return
    elif arguments.mode == "2":
        if arguments.initial_dirs:
            pedigree_tree_directory_paths = arguments.initial_dirs
            for directory in pedigree_tree_directory_paths:
                if not utility.verify_folder_path(directory):
                    print(f"Path {directory} does not exist or is not a directory")
                    return
            print(f"Running mode 2 selected with initial_dirs: {arguments.initial_dirs}")
            run_specific_error_pedigree_directories(paths=pedigree_tree_directory_paths,
                                                    error_pedigrees_folder_path=error_directory_path,
                                                    simulation_name=simulation_name)
        else:
            print("Mode 2 requires --initial-dirs.")
            return
    elif arguments.mode == "3":
        if arguments.initial_parent_dir:
            parent_paths_directory = arguments.initial_parent_dir
            if not utility.verify_folder_path(parent_paths_directory):
                print(f"Path {parent_paths_directory} does not exist or is not a directory")
                return
            pedigree_tree_directory_paths = get_absolute_paths_to_subfolders(parent_paths_directory)
            if not pedigree_tree_directory_paths:
                print(f"The directory {parent_paths_directory} does not contain any subdirectories")
                return
            print(f"Running mode 3 selected with initial_parent_dir: {parent_paths_directory}")
            run_specific_error_pedigree_directories(paths=pedigree_tree_directory_paths,
                                                    error_pedigrees_folder_path=error_directory_path,
                                                    simulation_name=simulation_name)
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
