from __future__ import annotations

import copy
import datetime
import os
import time
from builtins import int
from collections import defaultdict
from functools import partial
from pathlib import Path
from typing import Iterable

import igraph
import igraph as ig
from igraph import Flow
from lineagekit.core.CoalescentTree import CoalescentTree

from alignment.alignment_result import FullAlignmentResult, FailedClimbingAlignmentResult
from alignment.configuration import *
from alignment.potential_mrca_processed_graph import *
from alignment.subtree_matcher import *
from scripts.utility.alignment_utility import dict_has_duplicate_values


def get_subtree_matcher_for_coalescent_tree_proband(proband: int, proband_pedigree_ids: [int]):
    """
    Helper function returning the valid sub-alignments for a proband (which is simply one identity alignment).

    Args:
        proband (int): The id of the proband in the coalescent tree
        proband_pedigree_ids ([int]): The ids of the corresponding vertices in the pedigree
    """
    return {proband_pedigree_id: SubtreeMatcher(root_coalescent_tree=proband, root_pedigree=proband_pedigree_id)
            for proband_pedigree_id in proband_pedigree_ids}


def get_initial_simulation_mapping_for_mode(
        coalescent_tree: CoalescentTree,
        mode: ProbandInitialAssignmentsMode = default_proband_initial_assignments_mode):
    match mode:
        case ProbandInitialAssignmentsMode.PLOID:
            return {x: [x] for x in coalescent_tree.get_sink_vertices()}
        case ProbandInitialAssignmentsMode.INDIVIDUAL:
            return {x: [2 * (x // 2), (2 * (x // 2) + 1)] for x in coalescent_tree.get_sink_vertices()}
        case _:
            raise ValueError(f"Unknown initial matching mode: {mode}")


def get_pedigree_simulation_probands_for_alignment_mode(vertices: Iterable[int],
                                                        alignment_mode: ProbandInitialAssignmentsMode =
                                                        default_proband_initial_assignments_mode):
    # Returns the corresponding probands ids for a pedigree for the current matching mode
    match alignment_mode:
        case ProbandInitialAssignmentsMode.PLOID:
            return vertices
        case ProbandInitialAssignmentsMode.INDIVIDUAL:
            proband_individuals = {proband // 2 for proband in vertices}
            probands_all_ploids = [
                ploid
                for proband_individual in proband_individuals
                for ploid in [2 * proband_individual, 2 * proband_individual + 1]
            ]
            return probands_all_ploids
        case _:
            raise ValueError(f"Unknown alignment mode: {alignment_mode}")


source_vertex_label = "s"
target_vertex_label = "t"


class GraphMatcher:
    """
    This class runs the divide-and-conquer alignment algorithm on the given preprocessed pedigree and the coalescent
    tree.
    """

    class MatcherLogger:

        section_delimiter = 100 * "#"
        subsection_delimiter = 100 * '-'

        def __init__(self, logs_directory_path=logs_default_directory_name):
            """
            Initializes the logger and creates the log file.

            Args:
                logs_directory_path: The directory where the logs will be stored.
            """
            time_in_milliseconds = datetime.datetime.now().timestamp() * 1000
            log_filename = Path(logs_directory_path) / f"log_{int(time_in_milliseconds)}.txt"
            if not os.path.exists(logs_directory_path):
                os.makedirs(logs_directory_path)
            print(f"The log file path: {log_filename}")
            self.file = open(log_filename, 'w')

        def log(self, line: str):
            """
            Logs the passed line and flushes the buffer.

            Args:
                line (str): The line to log.
            """
            self.file.write(f"{line}\n")
            self.file.flush()

        def log_vertex_inference_time(self, coalescent_tree: CoalescentTree,
                                      spent_time: float, focal_vertex_coalescent_id: int, focal_vertex_children: [int],
                                      child_candidates_map: {int: [int]}, inference_result):
            """
            Logs the time spent on the coalescent vertex inference as well as the inference results themselves.

            Args:

                coalescent_tree (CoalescentTree): The coalescent tree.
                spent_time (float): The time spent on the inference.
                focal_vertex_coalescent_id (int): The coalescent vertex id for which the inference has been done.
                focal_vertex_children (list[int]): The coalescent vertex children vertices.
                child_candidates_map (dict[int, [int]]): Dictionary mapping every focal's child vertex to the list
                 of its pedigree candidates.
                inference_result: The list of candidates for the focal vertex.
            """
            self.log(f"Vertex level: {coalescent_tree.get_vertex_level(focal_vertex_coalescent_id)}")
            self.log(f"Solving the inference for {focal_vertex_coalescent_id} which has {len(focal_vertex_children)} "
                     f"children")
            self.log(f"The children's domain:")
            for child in focal_vertex_children:
                self.log(f"{child}: {len(child_candidates_map[child])}")
            self.log(f"Time taken for vertex inference: {spent_time}")
            total_candidates_found = len(inference_result)
            self.log(f"Total candidates found: {total_candidates_found}")

        def close(self):
            """
            Closes the log file.
            """
            self.file.close()

    def __init__(self, processed_graph: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree,
                 initial_mapping: dict[int, [int]], result_callback_function,
                 logs_path: str | Path = None,
                 alignment_vertex_mode: AlignmentVertexMode = default_alignment_vertex_mode,
                 alignment_edge_mode: AlignmentEdgeMode = default_alignment_edge_mode,
                 calculate_posterior_probabilities: bool = default_calculate_posterior_probabilities
                 ):
        """
        Initializes the GraphMatcher object.

        Args:
            processed_graph (PotentialMrcaProcessedGraph): The preprocessed pedigree graph that stores the
                "access vertices" through which a descendant vertex can reach the ancestor vertex.
            coalescent_tree (CoalescentTree): The coalescent tree with which the processed_graph is to be aligned.
            initial_mapping (dict[int, [int]]): The initial mapping between the proband vertices in the processed
                pedigree to the vertices in the coalescent tree.
            logs_path (str | Path): The path where the log file should be stored. None specifies that no logs
                should be stored.
            alignment_vertex_mode (AlignmentVertexMode): Specifies how extensive the search must be.
            alignment_edge_mode (AlignmentEdgeMode): Specifies whether edge alignments should be calculated.
        """
        self.pedigree = processed_graph
        self.coalescent_tree = coalescent_tree
        self.initial_mapping = initial_mapping
        self.alignment_vertex_mode = alignment_vertex_mode
        self.alignment_edge_mode = alignment_edge_mode
        self.result_callback_function = result_callback_function
        self.calculate_posterior_likelihoods = calculate_posterior_probabilities
        if logs_path and logs_enabled:
            self.logger = GraphMatcher.MatcherLogger(logs_directory_path=logs_path)
        else:
            self.logger = None

    def log_section_delimiter(self):
        self.log(GraphMatcher.MatcherLogger.section_delimiter)

    def log_subsection_delimiter(self):
        self.log(GraphMatcher.MatcherLogger.subsection_delimiter)

    def log(self, string: str):
        """
        Logs the given string if the logs are enabled. Does nothing otherwise.

        Args:
            string (str): The string to log.
        """
        if self.logger:
            self.logger.log(string)
        if print_enabled:
            print(string)

    def _find_potential_alignments(self):
        """
        This method finds the superset of all the valid alignments by climbing the tree and resolving all the
        'local' collisions. That is, all of these alignments won't have collisions on the parent-child level

        Returns:
             The set of the potential alignments.
        """
        coalescent_tree_vertex_to_subtrees: dict[int, dict[int, SubtreeMatcher]] = dict()
        clades = self.coalescent_tree.get_connected_components()
        if not clades:
            raise ValueError("The passed tree is empty")
        for clade in clades:
            max_level = max(self.coalescent_tree.get_vertex_level(v) for v in clade)
            clade_by_levels = [[] for _ in range(max_level + 1)]
            for v in clade:
                level = self.coalescent_tree.get_vertex_level(v)
                clade_by_levels[level].append(v)
            assert len(clade_by_levels[max_level]) == 1
            clade_root = clade_by_levels[max_level][0]
            for vertex in clade_by_levels[0]:
                coalescent_tree_vertex_to_subtrees[vertex] = get_subtree_matcher_for_coalescent_tree_proband(
                    vertex, self.initial_mapping[vertex])
            for level in clade_by_levels[1:]:
                proceed_climbing = True
                for vertex in level:
                    subtree_alignments = self._get_subtrees_from_children(
                        vertex,
                        coalescent_tree_vertex_to_subtrees
                    )
                    coalescent_tree_vertex_to_subtrees[vertex] = subtree_alignments
                    if len(subtree_alignments) == 0:
                        # Call the callback function with a failure alignment result object
                        focal_vertex_children: [int] = self.coalescent_tree.get_children(vertex)
                        child_to_pedigree_candidates = {
                            x: tuple(coalescent_tree_vertex_to_subtrees[x].keys())
                            for x in focal_vertex_children
                        }
                        failed_alignment_result = FailedClimbingAlignmentResult(
                            clade_root=clade_root,
                            failed_vertex=vertex,
                            child_to_pedigree_candidates=child_to_pedigree_candidates
                        )
                        # Notify the client code that climbing for the specified clade has failed
                        self.result_callback_function(failed_alignment_result)
                        proceed_climbing = False
                        break
                # Proceed to the next clade if climbing failed
                if not proceed_climbing:
                    break
        # tree_levels = self.coalescent_tree.get_levels()
        # if not tree_levels:
        #     raise ValueError("The passed tree is empty")
        # for vertex in tree_levels[0]:
        #     coalescent_tree_vertex_to_subtrees[vertex] = get_subtree_matcher_for_coalescent_tree_proband(
        #         vertex, self.initial_mapping[vertex])
        # for level in tree_levels[1:]:
        #     for vertex in level:
        #         subtree_alignments = self._get_subtrees_from_children(
        #             vertex,
        #             coalescent_tree_vertex_to_subtrees
        #         )
        #         if len(subtree_alignments) == 0:
        #             focal_vertex_children: [int] = self.coalescent_tree.get_children(vertex)
        #             child_to_pedigree_candidates = {
        #                 x: coalescent_tree_vertex_to_subtrees[x]
        #                 for x in focal_vertex_children
        #             }
        #             return FailClimbingResult(
        #                 tree_vertex_failure=vertex,
        #                 child_vertex_to_candidates=child_to_pedigree_candidates
        #             )
        #         coalescent_tree_vertex_to_subtrees[vertex] = subtree_alignments
        # Filtering the map, so that only the assignments for the roots of the clades remain
        coalescent_tree_vertex_to_subtrees = {
            x: value
            for x, value in coalescent_tree_vertex_to_subtrees.items()
            if x in self.coalescent_tree.get_founders()
        }
        # return SuccessClimbingResult(
        #     tree_vertex_to_subtree_alignments=coalescent_tree_vertex_to_subtrees
        # )
        return coalescent_tree_vertex_to_subtrees

    def find_alignments(self):
        """
        This method finds all the valid mapping between the given processed pedigree and
        every sub-clade within the coalescent tree. Every valid vertex alignment is passed to the callback function
        provided in the constructor
        """
        coalescent_tree_vertex_to_subtrees = self._find_potential_alignments()
        start_time = time.time()
        self.log_section_delimiter()
        self.log(f"Filtering the potential alignments. The alignment mode is {self.alignment_vertex_mode}")
        self._filter_alignments(coalescent_tree_vertex_to_subtrees)
        end_time = time.time()
        self.log(f"Time spent on building and filtering the results: {end_time - start_time}")
        # climb_result: ClimbingResult = self._find_potential_alignments()
        # match climb_result:
        #     case SuccessClimbingResult(tree_vertex_to_subtree_alignments=coalescent_tree_vertex_to_subtrees):
        #         start_time = time.time()
        #         self.log_section_delimiter()
        #         self.log(f"Filtering the potential alignments. The alignment mode is {self.alignment_vertex_mode}")
        #         self._filter_alignments(coalescent_tree_vertex_to_subtrees)
        #         end_time = time.time()
        #         self.log(f"Time spent on building and filtering the results: {end_time - start_time}")
        #     case FailClimbingResult(tree_vertex_failure=tree_vertex_failure,
        #                             child_vertex_to_candidates=child_vertex_to_candidates):
        #         pass

    def _verify_symmetry_holds(self, vertex_individual_id: int):
        """
        This function checks whether the ploids of the individual have the same children

        Args:
            vertex_individual_id: The individual id.

        Returns:
             Whether the individual's ploids have the same children.
        """
        paternal_ploid = 2 * vertex_individual_id
        maternal_ploid = paternal_ploid + 1
        return set(self.pedigree.get_children(paternal_ploid)) == set(self.pedigree.get_children(maternal_ploid))

    @staticmethod
    def _update_edge_alignment_for_symmetric_ploid(edge_alignment: dict[(int, int), [int]],
                                                   new_ploid: int, root_id: int,
                                                   coalescent_tree: CoalescentTree):
        other_ploid = PotentialMrcaProcessedGraph.get_other_ploid(new_ploid)
        children = coalescent_tree.get_children(root_id)
        for child in children:
            edge = (child, root_id)
            path: list[int] = edge_alignment[edge]
            removed_vertex = path.pop()
            assert removed_vertex == other_ploid
            path.append(new_ploid)

    def _get_symmetric_alignment_result(self, alignment_result: FullAlignmentResult):
        # Using the results obtained for the other ploid
        symmetric_result = copy.deepcopy(alignment_result)
        clade_root = symmetric_result.clade_root
        old_ploid = symmetric_result[clade_root]
        new_ploid = self.pedigree.get_other_ploid(old_ploid)
        # Update the vertex alignment
        symmetric_result.vertex_alignment[clade_root] = new_ploid
        edge_alignments_to_update = []
        if symmetric_result.example_edge_alignment:
            edge_alignments_to_update.append(symmetric_result.example_edge_alignment)
        elif symmetric_result.edge_alignments:
            edge_alignments_to_update = symmetric_result.edge_alignments
        for edge_alignment in edge_alignments_to_update:
            self._update_edge_alignment_for_symmetric_ploid(edge_alignment=edge_alignment,
                                                            root_id=clade_root,
                                                            new_ploid=new_ploid,
                                                            coalescent_tree=self.coalescent_tree
                                                            )
        # Update the posterior probabilities if calculated
        if symmetric_result.posterior_probabilities:
            symmetric_result.posterior_probabilities.vertex_posterior_probabilities_for_vertex_alignment.pop(old_ploid)
            symmetric_result.posterior_probabilities.vertex_posterior_probabilities_for_vertex_alignment[new_ploid] = 1
        return symmetric_result

    def _filter_alignments(self, coalescent_tree_vertex_to_subtrees: dict):
        match self.alignment_vertex_mode:
            case AlignmentVertexMode.ALL_ALIGNMENTS:
                root_matcher_filter_function = self._get_verified_alignments_brute_force
            case AlignmentVertexMode.EXAMPLE_PER_ROOT_ASSIGNMENT:
                root_matcher_filter_function = self._find_example_alignment_for_root_assignment
            case _:
                raise ValueError("Unsupported AlignmentExtentMode")
        original_callback = self.result_callback_function
        # Iterating over clades
        for clade_root_vertex, clade_root_vertex_assignments in coalescent_tree_vertex_to_subtrees.items():
            root_vertex_matchers: [SubtreeMatcher] = clade_root_vertex_assignments.values()
            log_message = (f"Filtering the alignments for the root {clade_root_vertex}\n"
                           f"There are {len(root_vertex_matchers)} potential assignments for this root")
            self.log(log_message)
            pedigree_individuals_respecting_ploid_symmetry = set()
            root_assignments_to_process = []
            # Since the two ploids of the same individual are symmetrical
            # (unless the pedigree has been modified after parsing), we can only process one ploid and use
            # its results for the other
            # This check is performed because the user can modify the pedigree in such a way that the symmetry doesn't
            # hold
            for root_matcher in root_vertex_matchers:
                pedigree_candidate = root_matcher.root_pedigree
                pedigree_individual = pedigree_candidate // 2
                if pedigree_individual in pedigree_individuals_respecting_ploid_symmetry:
                    self.log(f"Using the symmetry of the pedigree for the pedigree candidate {pedigree_candidate}")
                    continue
                root_assignments_to_process.append(root_matcher)
                if self._verify_symmetry_holds(vertex_individual_id=pedigree_individual):
                    pedigree_individuals_respecting_ploid_symmetry.add(pedigree_individual)

            def intermediary_callback(alignment_result: FullAlignmentResult):
                # This callback takes advantage of the symmetry in the pedigree to speed up the calculations
                root_pedigree_candidate = alignment_result[alignment_result.clade_root]
                root_pedigree_candidate_individual = root_pedigree_candidate // 2
                if root_pedigree_candidate_individual in pedigree_individuals_respecting_ploid_symmetry:
                    symmetric_result = self._get_symmetric_alignment_result(alignment_result)
                    original_callback(symmetric_result)
                original_callback(alignment_result)

            self.result_callback_function = intermediary_callback
            # Processing all the selected pedigree candidates for the clade root
            time_start_total = time.time()
            for root_matcher in root_assignments_to_process:
                self.log_subsection_delimiter()
                self.log(f"Filtering the alignments for the pedigree candidate {root_matcher.root_pedigree}")
                validation_start = time.time()
                root_matcher_filter_function(root_matcher)
                validation_end = time.time()
                time_taken = validation_end - validation_start
                self.log(f"Time spent on this root candidate: {time_taken}")
            time_end_total = time.time()
            time_taken_total = time_end_total - time_start_total
            self.log(f"Total time spent on filtering: {time_taken_total}")
        self.result_callback_function = original_callback

    def _get_verified_alignments_brute_force(self, matcher: SubtreeMatcher):
        alignments = matcher.get_all_subtree_alignments()
        self.log(f"There are {len(alignments)} alignments for the matcher, filtering the results")
        valid_alignments = 0
        for alignment in alignments:
            alignment_verification_result = self._verify_valid_alignment(potential_alignment=alignment)
            if alignment_verification_result:
                self.result_callback_function(alignment_verification_result)
                valid_alignments += 1
        self.log(f"The number of valid alignments: {valid_alignments}")

    def _find_example_alignment_for_root_assignment(self, root_matcher: SubtreeMatcher):
        self.log(f"Finding an example alignment for the assignment "
                 f"{root_matcher.root_coalescent_tree}: {root_matcher.root_pedigree}")
        alignments = root_matcher.get_all_subtree_alignments()
        self.log(f"There are {len(alignments)} potential alignments")
        for alignment in alignments:
            alignment_verification = self._verify_valid_alignment(potential_alignment=alignment)
            if alignment_verification:
                self.result_callback_function(alignment_verification)
                break

    @staticmethod
    def _are_paths_disjoint(paths):
        seen = set()
        for path in paths:
            for vertex in path:
                if vertex in seen:
                    return False
                seen.add(vertex)
        return True

    @staticmethod
    def _map_trimmed_path_to_readable_format(vertex_alignment, edge, trimmed_path, network_graph):
        _, parent = edge
        readable_path = []
        for index, vertex in enumerate(trimmed_path):
            if index % 2 == 1:
                # Every second vertex is an access vertex which shouldn't be reported
                continue
            pedigree_ploid_id = int(network_graph.vs[vertex]["name"])
            readable_path.append(pedigree_ploid_id)
        last_path_vertex = vertex_alignment[parent]
        readable_path.append(last_path_vertex)
        return readable_path

    @staticmethod
    def _map_edge_alignment_to_meaningful_ids(edge_alignment, vertex_alignment, network_graph):
        readable_edge_alignment = dict()
        for edge, path in edge_alignment:
            readable_path = GraphMatcher._map_trimmed_path_to_readable_format(vertex_alignment, edge,
                                                                              path, network_graph)
            readable_edge_alignment[edge] = readable_path
        return readable_edge_alignment

    def _edge_alignment_search_simple_product(self, edge_to_all_pedigree_paths: dict):
        edge_to_path_tuples = [
            [(edge, path) for path in paths]
            for edge, paths in edge_to_all_pedigree_paths.items()
        ]
        # This is a generator (lazy evaluation)
        possible_edge_alignments = product(*edge_to_path_tuples)
        valid_edge_alignments = []
        for possible_edge_alignment in possible_edge_alignments:
            paths = [path for _, path in possible_edge_alignment]
            if self._are_paths_disjoint(paths):
                valid_edge_alignments.append(possible_edge_alignment)
                if self.alignment_vertex_mode == AlignmentVertexMode.EXAMPLE_PER_ROOT_ASSIGNMENT:
                    break
        return valid_edge_alignments

    def _edge_alignment_search_sort_by_value_space(self, edge_to_all_pedigree_paths: dict):
        # Sort edges by number of paths (most constrained first)
        sorted_edges = sorted(
            edge_to_all_pedigree_paths.items(),
            key=lambda item: len(item[1])
        )
        edge_number = len(sorted_edges)
        valid_edge_alignments = []

        def validate_valid_edge_alignment_extension(edge_alignment, new_path):
            # Avoid constructing the set if the edge alignment is empty
            if not edge_alignment:
                return True
            new_path_vertices = frozenset(new_path)
            for _, path in edge_alignment:
                if new_path_vertices.intersection(path):
                    return False
            return True

        def backtrack_search(edge_index: int, current_edge_alignment: list):
            if edge_index == edge_number:
                valid_edge_alignments.append(current_edge_alignment.copy())
                return
            edge, paths = sorted_edges[edge_index]
            for path in paths:
                if validate_valid_edge_alignment_extension(current_edge_alignment, path):
                    current_edge_alignment.append((edge, path))
                    backtrack_search(edge_index + 1, current_edge_alignment)
                    current_edge_alignment.pop()
                    if (len(valid_edge_alignments) > 1 and
                            self.alignment_vertex_mode == AlignmentVertexMode.EXAMPLE_PER_ROOT_ASSIGNMENT):
                        return

        backtrack_search(0, [])
        return valid_edge_alignments

    def _build_edge_to_path_map(self, alignment: dict, network_graph: igraph.Graph):
        edge_to_all_pedigree_paths = dict()
        for parent, child in self.coalescent_tree.edges:
            child_pedigree = alignment.get(child, None)
            parent_pedigree = alignment.get(parent, None)
            if child_pedigree is None or parent_pedigree is None:
                continue
            parent_pedigree_access_vertex = f"{parent_pedigree}'"
            child_pedigree_id = network_graph.vs.find(name=str(child_pedigree)).index
            parent_pedigree_id = network_graph.vs.find(name=str(parent_pedigree)).index
            parent_pedigree_access_vertex_id = network_graph.vs.find(name=parent_pedigree_access_vertex).index
            simple_paths: list[list[int]] = network_graph.get_all_simple_paths(v=child_pedigree_id,
                                                                               to=parent_pedigree_id,
                                                                               mode=igraph.OUT)
            for simple_path in simple_paths:
                simple_path.remove(parent_pedigree_id)
                # We can remove the access vertex as well!
                # Notice that any vertex that goes through the access vertex must go through the original
                # vertex as well
                # Now, consider the case where the parent vertex is not the candidate for the root vertex. In this case,
                # it will be present in the path going from the parent vertex to its parent. Therefore, if any other
                # vertex goes through the parent vertex's access vertex, we will have a collision
                # On the other hand, if the parent vertex is the candidate for the root of the clade, then no other
                # vertex can actually reach it (as this would mean that there is a cycle in the pedigree)
                simple_path.remove(parent_pedigree_access_vertex_id)
            edge_to_all_pedigree_paths[(child, parent)] = simple_paths
        return edge_to_all_pedigree_paths

    def _get_possible_edge_alignment_generator(self, alignment: dict, network_graph: igraph.Graph):
        edge_to_all_pedigree_paths = self._build_edge_to_path_map(alignment, network_graph)
        edge_to_path_tuples = [
            [(edge, path) for path in paths]
            for edge, paths in edge_to_all_pedigree_paths.items()
        ]
        # This is a generator (lazy evaluation)
        return product(*edge_to_path_tuples)

    def _get_edge_alignments_for_vertex_alignment(self, potential_alignment: FullAlignmentResult,
                                                  network_graph: igraph.Graph,
                                                  search_function) -> FullAlignmentResult:
        alignment = potential_alignment.vertex_alignment
        edge_to_all_pedigree_paths = self._build_edge_to_path_map(alignment=alignment, network_graph=network_graph)
        edge_alignments = search_function(edge_to_all_pedigree_paths)
        partial_func = partial(self._map_edge_alignment_to_meaningful_ids,
                               vertex_alignment=alignment,
                               network_graph=network_graph)
        valid_edge_alignments = list(map(partial_func, edge_alignments))
        if len(valid_edge_alignments) > 0:
            potential_alignment.edge_alignments = valid_edge_alignments
            potential_alignment.is_valid = True
            if self.calculate_posterior_likelihoods:
                potential_alignment.calculate_vertex_inclusion_probabilities()
        return potential_alignment

    def _get_example_edge_alignment(self, potential_alignment: FullAlignmentResult,
                                    network_graph: igraph.Graph) -> FullAlignmentResult:
        alignment = potential_alignment.vertex_alignment
        potential_edge_alignments = self._get_possible_edge_alignment_generator(alignment=alignment,
                                                                                network_graph=network_graph)
        for possible_edge_alignment in potential_edge_alignments:
            paths = [path for _, path in possible_edge_alignment]
            if self._are_paths_disjoint(paths):
                potential_alignment.example_edge_alignment = possible_edge_alignment
                return potential_alignment

    def _verify_example_edge_alignment_from_flow(self, alignment: dict[int, int],
                                                 alignment_flow: Flow) -> dict[(int, int), [int]]:
        """
        Creates an example edge alignment from the obtained flow

        Args:
            alignment: The vertex alignment
            alignment_flow: The corresponding maximum flow

        Returns:
             The edge alignment if the flow produces a correct alignment or None otherwise
        """

        def get_edge_flow(edge_child: int, edge_parent: int):
            edge_id = graph.get_eid(edge_child, edge_parent)
            return alignment_flow.flow[edge_id]

        def get_pedigree_path(descendant: int, ancestor: int):
            current_path = [descendant]
            current_vertex = descendant
            current_vertex_id = graph.vs.find(name=str(current_vertex)).index
            ancestor_id = graph.vs.find(name=str(ancestor)).index
            target_id = graph.vs.find(name=target_vertex_label).index
            while current_vertex_id != ancestor_id:
                neighbours = graph.neighbors(vertex=current_vertex_id, mode="out")
                neighbours = [x for x in neighbours if get_edge_flow(edge_child=current_vertex_id, edge_parent=x) > 0]
                try:
                    neighbours.remove(target_id)
                except ValueError:
                    pass
                if not neighbours:
                    # The flow doesn't connect the descendant with the ancestor, the flow doesn't produce the alignment
                    raise ValueError("The flow doesn't produce a valid alignment")
                assert len(neighbours) == 1
                # The vertex that we've reached is the access vertex of the next pedigree individual on the path
                access_vertex_id = neighbours[0]
                access_vertex_neighbors = graph.neighbors(vertex=access_vertex_id, mode="out")
                # The access vertex is only connected to its corresponding original vertex
                assert len(access_vertex_neighbors) == 1
                next_vertex_id = access_vertex_neighbors[0]
                next_vertex_name = int(graph.vs[next_vertex_id]["name"])
                current_path.append(next_vertex_name)
                current_vertex_id = next_vertex_id
            return current_path

        graph: ig.Graph = alignment_flow.graph
        tree_non_leaf_levels = self.coalescent_tree.get_levels()[1:]
        edge_alignment = dict()
        for tree_level in tree_non_leaf_levels:
            for vertex in tree_level:
                vertex_pedigree = alignment.get(vertex, None)
                # If the vertex is not in the alignment (for example, it's a part of another clade, skip it)
                if not vertex_pedigree:
                    continue
                vertex_children = self.coalescent_tree.get_children(vertex)
                for child in vertex_children:
                    child_pedigree = alignment.get(child, None)
                    if not child_pedigree:
                        continue
                    try:
                        pedigree_path = get_pedigree_path(descendant=child_pedigree, ancestor=vertex_pedigree)
                    except ValueError:
                        return None
                    edge_alignment[(child, vertex)] = pedigree_path
        return edge_alignment

    def _run_maximum_flow(self, potential_alignment: FullAlignmentResult, network_graph: igraph.Graph,
                          proband_number: int) -> FullAlignmentResult:
        # Verify that the default value is False
        assert not potential_alignment.is_valid
        match self.alignment_edge_mode:
            case AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT:
                maximum_flow = network_graph.maxflow(source=source_vertex_label, target=target_vertex_label,
                                                     capacity='capacity')
                maximum_flow_value = maximum_flow.value
                necessary_condition = maximum_flow_value == proband_number
                # If the necessary condition doesn't hold, there are no edge alignments
                if not necessary_condition:
                    return potential_alignment
                # We can additionally check if the alignment is correct
                # If so, we can use it as an example alignment
                edge_alignment = self._verify_example_edge_alignment_from_flow(
                    alignment=potential_alignment.vertex_alignment,
                    alignment_flow=maximum_flow
                )
                if edge_alignment is not None:
                    potential_alignment.example_edge_alignment = edge_alignment
                    potential_alignment.is_valid = True
                    return potential_alignment
                # Otherwise, we need to perform a search until we find a solution or run out of options
                return self._get_edge_alignments_for_vertex_alignment(
                    network_graph=network_graph,
                    potential_alignment=potential_alignment,
                    search_function=self._edge_alignment_search_simple_product
                )
            case AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS:
                maximum_flow_value = network_graph.maxflow_value(source=source_vertex_label,
                                                                 target=target_vertex_label,
                                                                 capacity='capacity')
                is_valid = maximum_flow_value == proband_number
                if not is_valid:
                    return potential_alignment
                return self._get_edge_alignments_for_vertex_alignment(
                    network_graph=network_graph,
                    potential_alignment=potential_alignment,
                    search_function=self._edge_alignment_search_simple_product
                )

    def _verify_valid_alignment(self, potential_alignment: dict[int, int]) -> FullAlignmentResult:
        """
        Verifies that the given (possibly partial) alignment is correct.

        Args:
            potential_alignment (dict): The alignment to be verified.

        Returns:
            The alignment result object calculated based on the specified problem mode.
        """
        if dict_has_duplicate_values(potential_alignment):
            return FullAlignmentResult(
                vertex_alignment=potential_alignment,
                is_valid=False,
                clade_root=None
            )
        clade_root = None
        network_graph = ig.Graph(directed=True)
        vertices_to_add = set()
        edges_to_add = set()
        edge_capacities = dict()
        vertices_to_add.update([source_vertex_label, target_vertex_label])
        proband_number = 0
        for tree_vertex in potential_alignment:
            tree_vertex_pedigree = potential_alignment[tree_vertex]
            tree_vertex_pedigree_str = str(tree_vertex_pedigree)
            vertices_to_add.add(tree_vertex_pedigree_str)
            tree_vertex_parents = [x for x in self.coalescent_tree.get_parents(tree_vertex) if x in potential_alignment]
            # If the vertex is a proband, connect it to the source
            children = self.coalescent_tree.get_children(tree_vertex)
            if not children:
                edges_to_add.add((source_vertex_label, tree_vertex_pedigree_str))
                proband_number += 1
                if not tree_vertex_parents:
                    assert clade_root is None
                    clade_root = tree_vertex
                    edges_to_add.add((tree_vertex_pedigree_str, target_vertex_label))
                continue
            assert len(children) > 0
            children_pedigree = [potential_alignment[x] for x in children if potential_alignment.get(x) is not None]
            # The number of children present in the alignment
            children_number = len(children_pedigree)
            vertices, edges, capacities = self._add_edges_to_mrca_from_descendants(tree_vertex_pedigree,
                                                                                   children_pedigree)

            vertices_to_add.update(vertices)
            edges_to_add.update(edges)
            assert not [x for x in capacities if x in edge_capacities]
            edge_capacities.update(capacities)
            edge = tree_vertex_pedigree_str, target_vertex_label
            edges_to_add.add(edge)
            if not tree_vertex_parents:
                assert clade_root is None
                clade_root = tree_vertex
                edge_capacities[edge] = children_number
            else:
                edge_capacities[edge] = children_number - 1
        vertices_to_add = list(vertices_to_add)
        network_graph.add_vertices(vertices_to_add)
        edges = list(edges_to_add)
        capacities = [edge_capacities.get(edge, 1) for edge in edges]
        network_graph.add_edges(edges, attributes={"capacity": capacities})
        assert clade_root is not None
        alignment_result = FullAlignmentResult(
            vertex_alignment=potential_alignment,
            clade_root=clade_root,
            is_valid=False
        )
        return self._run_maximum_flow(potential_alignment=alignment_result, network_graph=network_graph,
                                      proband_number=proband_number)

    def _add_edges_to_mrca_from_descendants(self, mrca: int, descendant_vertices: [int]):
        mrca_access_label = f"{mrca}'"
        vertices_to_add = set()
        edges_to_add = set()
        edge_to_capacity = dict()
        vertices_to_add.add(mrca_access_label)
        edge_to_add = (mrca_access_label, str(mrca))
        edges_to_add.add(edge_to_add)
        edge_to_capacity[edge_to_add] = len(descendant_vertices)
        for descendant in descendant_vertices:
            current_level_vertices = [mrca]
            while current_level_vertices:
                next_level_vertices = set()
                for vertex in current_level_vertices:
                    if vertex == descendant:
                        continue
                    neighbors = self.pedigree.vertex_to_ancestor_map[descendant].get(vertex, ())
                    next_level_vertices.update(neighbors)
                    vertex_str = str(vertex)
                    vertex_access_label = f"{vertex_str}'"
                    vertices_to_add.update((vertex_access_label, vertex_str))
                    neighbors = [str(neighbor) for neighbor in neighbors]
                    vertices_to_add.update(neighbors)
                    if vertex != mrca:
                        edges_to_add.add((vertex_access_label, vertex_str))
                    edges = [(neighbor, vertex_access_label) for neighbor in neighbors]
                    edges_to_add.update(edges)
                current_level_vertices = next_level_vertices
        return vertices_to_add, edges_to_add, edge_to_capacity

    def _get_subtrees_from_children(self, focal_vertex: int, vertex_subtree_dict: {int: {int: SubtreeMatcher}}):
        """
        Finds all the valid alignments for the sub-clade where the given focal_vertex is the root.

        Args:
            focal_vertex: The root of the sub-clade for which the inference is to be performed.
            vertex_subtree_dict (dict[int, [int]]): The dictionary containing all the valid alignment for
             the focal vertex.

        Returns:
             The list of all the valid sub-clade alignments for the specified focal vertex. The resulting
             list consists of SubtreeMatcher subtree_matcher.SubtreeMatcher objects.
        """
        focal_vertex_children: [int] = self.coalescent_tree.get_children(focal_vertex)
        if not focal_vertex_children:
            raise Exception("Isolated vertex in the coalescent tree")
        if len(focal_vertex_children) == 1:
            child, = focal_vertex_children
            result = {x: SubtreeMatcher(root_coalescent_tree=focal_vertex, root_pedigree=x,
                                        children_assignments=[{child: y}]) for x, y in
                      vertex_subtree_dict[child].items()}
            return result
        if print_enabled:
            print(section_separator)
            print(f"Inference for {focal_vertex} ({self.coalescent_tree.get_vertex_level(focal_vertex)}),"
                  f" there are {len(focal_vertex_children)} children")
            for child in focal_vertex_children:
                print(f"{child}: {len(vertex_subtree_dict[child])}")
        inference_start = time.time()
        (children_order, inference_result) = self._get_pmracs_for_vertices(
            coalescent_vertex_to_candidates=vertex_subtree_dict,
            vertices_coalescent_ids=focal_vertex_children
        )
        inference_end = time.time()
        time_taken = inference_end - inference_start
        self.log(f"Inference time: {time_taken}")
        self.log_subsection_delimiter()
        if self.logger:
            self.logger.log_vertex_inference_time(self.coalescent_tree, time_taken, focal_vertex, focal_vertex_children,
                                                  vertex_subtree_dict, inference_result)
        result_build_start = time.time()
        candidate_subtree_matcher_dictionary = dict()
        inference_result_keys = list(inference_result)
        for focal_vertex_candidate in inference_result_keys:
            children_assignments = inference_result[focal_vertex_candidate]
            children_dictionaries = []
            for children_assignment in children_assignments:
                assert len(children_order) == len(children_assignment) == len(focal_vertex_children)
                children_dictionary = dict()
                for index, child_assigned_pedigree_id in enumerate(children_assignment):
                    child_coalescent_id = children_order[index]
                    assigned_child_matcher = vertex_subtree_dict[child_coalescent_id][child_assigned_pedigree_id]
                    children_dictionary[child_coalescent_id] = assigned_child_matcher
                children_dictionaries.append(children_dictionary)
            candidate_subtree_matcher_dictionary[focal_vertex_candidate] = SubtreeMatcher(
                root_coalescent_tree=focal_vertex,
                root_pedigree=focal_vertex_candidate,
                children_assignments=children_dictionaries)
            inference_result.pop(focal_vertex_candidate)
        result_build_end = time.time()
        self.log(f"Building the result {result_build_end - result_build_start}")
        self.log(f"There are {len(candidate_subtree_matcher_dictionary)} resulting assignments")
        return candidate_subtree_matcher_dictionary

    # ----------------------------------------- Alignment logic ------------------------------------------------------

    def _verify_pmrca_for_vertex_pair(self, first_vertex: int, second_vertex: int, common_ancestor: int):
        """
        Verifies that the given vertex is a potential mrca for the given pair of vertices.

        Args:
            first_vertex (int): The first vertex from the pair.
            second_vertex (int): The second vertex from the pair.
            common_ancestor (int): A common ancestor of the given pair of vertices that is being verified.

        Returns:
            True if the given vertex is potential mrca for the given pair of vertices, False otherwise.
        """
        first_vertex_access_vertices = self.pedigree.vertex_to_ancestor_map[first_vertex][common_ancestor]
        return not (len(first_vertex_access_vertices) == 1 and
                    first_vertex_access_vertices ==
                    self.pedigree.vertex_to_ancestor_map[second_vertex][common_ancestor])

    def _triplet_condition_holds(self, first_vertex: int, second_vertex: int, third_vertex: int,
                                 potential_mrca_candidate: int):
        """
        Verifies that the potential MRCA candidate satisfies all the conditions involving the
        third vertex. More specifically, it verifies that the `Hall's condition
        <https://en.wikipedia.org/wiki/Hall%27s_marriage_theorem>`_ is satisfied for the
        triplet (first_vertex, second_vertex, third_vertex).

        Args:
            first_vertex: The first vertex from the triplet.
            second_vertex: The second vertex from the triplet.
            third_vertex: The third vertex from the triplet.
            potential_mrca_candidate: A common ancestor of the given triplet of vertices that satisfies the
             Hall's condition for the (first_vertex, second_vertex) pair.

        Returns:
             True if potential mrca candidate is indeed a potential mrca, False otherwise.
        """
        first_vertex_access = self.pedigree.vertex_to_ancestor_map[first_vertex][potential_mrca_candidate]
        if len(first_vertex_access) > 2:
            return True
        second_vertex_access = self.pedigree.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
        current_access = set(first_vertex_access)
        # Add elements from the second access and check the distinct count
        for vertex in second_vertex_access:
            current_access.add(vertex)
            if len(current_access) > 2:
                return True
        third_vertex_access = self.pedigree.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
        # Add elements from the third access and check the distinct count
        for vertex in third_vertex_access:
            current_access.add(vertex)
            if len(current_access) > 2:
                return True
        # If the loop completes, there are 2 or fewer distinct elements
        return False

    def _filter_common_ancestors_for_vertex_pair(self, first_candidate: int, second_candidate: int):
        """
        Finds the pmracs vertices for the given vertex pair. Important: This method does not check whether the passed
        vertices are different.

        Args:
            first_candidate (int): The first vertex id.
            second_candidate (int): The second vertex id.

        Returns:
             The list of pmracs for the given vertices.
        """
        verified_ancestors = []
        if not self.pedigree.has_parents(first_candidate) or not self.pedigree.has_parents(second_candidate):
            return verified_ancestors
        while self.pedigree.get_parents(first_candidate) == self.pedigree.get_parents(second_candidate):
            verified_ancestors.extend(self.pedigree.get_parents(first_candidate))
            [first_candidate, second_candidate] = self.pedigree.get_parents(first_candidate)
            if not self.pedigree.has_parents(first_candidate) or not self.pedigree.has_parents(second_candidate):
                return verified_ancestors
        first_ancestors = self.pedigree.vertex_to_ancestor_map[first_candidate]
        second_ancestors = self.pedigree.vertex_to_ancestor_map[second_candidate]
        if len(second_ancestors) < len(first_ancestors):
            first_ancestors, second_ancestors = second_ancestors, first_ancestors
        verified_ancestors.extend([ancestor for ancestor in first_ancestors if ancestor in second_ancestors and
                                   self._verify_pmrca_for_vertex_pair(first_candidate, second_candidate, ancestor)
                                   ])
        return verified_ancestors

    def _get_pmracs_for_vertex_pair(self, first: int, second: int, coalescent_vertex_to_candidates: {int: [int]}):
        """
        Explores the possible candidates for the given coalescent vertices and finds the corresponding pmracs.

        Args:
            first: The first coalescent vertex id.
            second: The second coalescent vertex id.
            coalescent_vertex_to_candidates (dict[int, [int]]): A dictionary mapping coalescent vertex ids to the list
             of pedigree vertices.

        Returns:
            The dictionary mapping a valid pmrca id to the corresponding pedigree assignments for the passed
            coalescent vertices.
        """
        result = defaultdict(list)
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        if len(first_vertex_candidates) > len(second_vertex_candidates):
            first_vertex_candidates, second_vertex_candidates = second_vertex_candidates, first_vertex_candidates
        for first_vertex_candidate in first_vertex_candidates:
            for second_vertex_candidate in second_vertex_candidates:
                if first_vertex_candidate == second_vertex_candidate:
                    continue
                verified_ancestors = self._filter_common_ancestors_for_vertex_pair(
                    first_vertex_candidate,
                    second_vertex_candidate)
                for verified_ancestor in verified_ancestors:
                    result[verified_ancestor].append((first_vertex_candidate, second_vertex_candidate))
        return result

    def _get_pmracs_for_vertex_triple_iterative(self, first: int, second: int, third: int,
                                                coalescent_vertex_to_candidates: {int: [int]}):
        """
        Explores the possible candidates for the given coalescent vertices and finds the corresponding pmracs.

        Args:
            first: The first coalescent vertex.
            second:The second coalescent vertex.
            third: The third coalescent vertex.
            coalescent_vertex_to_candidates (dict[int, [int]]): A dictionary mapping coalescent vertex ids to the list
             of pedigree vertices.

        Returns:
             The dictionary mapping a valid pmrca id to the corresponding pedigree assignments for the passed
             coalescent vertices.
        """
        triplet_cache = dict()

        def get_triplet_tuple(first_element, second_element, third_element):
            triplet_tuple = (first_element, second_element, third_element)
            if triplet_tuple in triplet_cache:
                return triplet_cache[triplet_tuple]
            triplet_cache[triplet_tuple] = triplet_tuple
            return triplet_tuple

        first_second_pair_result = self._get_pmracs_for_vertex_pair(first, second, coalescent_vertex_to_candidates)
        first_third_pair_result = self._get_pmracs_for_vertex_pair(first, third, coalescent_vertex_to_candidates)
        first_second_dict = defaultdict(dict)
        first_third_dict = defaultdict(dict)
        for (pair_result, dictionary) in ((first_second_pair_result, first_second_dict),
                                          (first_third_pair_result, first_third_dict)):
            for valid_assignment in pair_result.items():
                (pmrca_candidate, children_assignments) = valid_assignment
                for children_assignment in children_assignments:
                    (first_candidate, other_candidate) = children_assignment
                    if first_candidate not in dictionary[pmrca_candidate]:
                        dictionary[pmrca_candidate][first_candidate] = []
                    dictionary[pmrca_candidate][first_candidate].append(other_candidate)
        # Letting the GC free the memory if necessary
        del first_second_pair_result
        del first_third_pair_result
        shared_common_ancestors = first_second_dict.keys() & first_third_dict.keys()
        result = defaultdict(list)
        for shared_common_ancestor in shared_common_ancestors:
            shared_first_vertex_assignments = (first_second_dict[shared_common_ancestor].keys() &
                                               first_third_dict[shared_common_ancestor].keys())
            for first_candidate in shared_first_vertex_assignments:
                for second_candidate in first_second_dict[shared_common_ancestor][first_candidate]:
                    for third_candidate in first_third_dict[shared_common_ancestor][first_candidate]:
                        if second_candidate == third_candidate:
                            continue
                        # TODO: Consider adding the common-parents speed-up
                        current_triplet_tuple = get_triplet_tuple(first_candidate, second_candidate, third_candidate)
                        if (self._verify_pmrca_for_vertex_pair(second_candidate, third_candidate,
                                                               shared_common_ancestor)
                                and self._triplet_condition_holds(first_candidate, second_candidate,
                                                                  third_candidate, shared_common_ancestor)):
                            result[shared_common_ancestor].append(current_triplet_tuple)
        return result

    def _get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
                                 coalescent_vertex_to_candidates: {int: [int]}):
        """
        This function calculates the potential most recent common ancestor (pmrca) for the given vertices.

        Args:
            vertices_coalescent_ids ([int]): The ids for the coalescent vertices.
            coalescent_vertex_to_candidates (dict[int, [int]]): A dictionary mapping an id of a coalescent vertex to the
             list of ids of pedigree vertices that can be assigned to this coalescent vertex.

        Returns:
             All the valid pmracs for the given vertices.
        """
        vertices_length = len(vertices_coalescent_ids)
        vertices_coalescent_ids = sorted(vertices_coalescent_ids,
                                         key=lambda child: len(coalescent_vertex_to_candidates[child]),
                                         reverse=False
                                         )
        if vertices_length == 2:
            result = self._get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                      coalescent_vertex_to_candidates)
        elif vertices_length == 3:
            result = self._get_pmracs_for_vertex_triple_iterative(vertices_coalescent_ids[0],
                                                                  vertices_coalescent_ids[1],
                                                                  vertices_coalescent_ids[2],
                                                                  coalescent_vertex_to_candidates)
        else:
            result = self._get_pmracs_for_vertices_dfs(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        return vertices_coalescent_ids, result

    def _get_pmracs_for_vertices_dfs(self, vertices_coalescent_ids: [int],
                                     coalescent_vertex_to_candidates: {int: [int]}):
        """
        Explores the search graph in a dfs manner storing the results only for the currently explored branch.

        Args:
            vertices_coalescent_ids (list[int]): The ids for the coalescent vertices.
            coalescent_vertex_to_candidates (dict[int, [int]]): A dictionary mapping an id of a coalescent vertex to
             the list of ids of pedigree vertices that can be assigned to this coalescent vertex.

        Returns:
            All the valid pmracs for the given vertices.
        """
        result = []
        vertices_length = len(vertices_coalescent_ids)
        # Specifies the index of the currently explored coalescent vertex
        current_index = 0
        # Specifies the currently explored partial assignment
        current_partial_assignment = dict()
        partial_result = dict()
        # The current_stack variable is a stack keeping track of the pedigree ploid assignments that are to be
        # explored for the current coalescent vertex
        # The stack has a format like this:
        # {pedigree candidates for the first vertex} -1 {pedigree candidates for the second vertex} -1 ...
        current_stack = []
        # Adding the pedigree candidates for the first tree vertex
        current_stack.extend(coalescent_vertex_to_candidates[vertices_coalescent_ids[0]])
        # A cache used for reducing the redundant creation of identical sets
        set_cache = dict()

        def cache_set(set_to_cache: frozenset):
            if set_to_cache in set_cache:
                return set_cache[set_to_cache]
            set_cache[set_to_cache] = set_to_cache
            return set_to_cache

        def free_current_vertex_data():
            current_coalescent_vertex = vertices_coalescent_ids[current_index]
            current_fixed_vertex = current_partial_assignment[current_coalescent_vertex]
            partial_result.pop(current_fixed_vertex)
            current_partial_assignment.pop(current_coalescent_vertex)

        while current_stack:
            vertex = current_stack.pop()
            if vertex in partial_result:
                continue
            # This signals that we've explored all the possible assignments for the current coalescent vertex
            if vertex == -1:
                current_index -= 1
                free_current_vertex_data()
                continue
            vertex_ancestors = self.pedigree.get_vertex_ancestors(vertex)
            verified_assignments = None
            if not current_index:
                # We have only mapped the first tree vertex
                first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(vertex)
                                             if len(self.pedigree.get_children(x)) >= vertices_length]
                if first_candidate_ancestors:
                    verified_assignments = (
                        (vertex,),
                        ([(ancestor, {cache_set(frozenset(self.pedigree.vertex_to_ancestor_map[vertex][ancestor])): 1})
                          for ancestor in vertex_ancestors])
                    )
            else:
                # Applying the new constraints for the next tree vertex using the Hall's theorem
                new_verified_candidates_list = []
                previously_assigned_vertex = current_partial_assignment[vertices_coalescent_ids[current_index - 1]]
                (assigned_children, verified_candidates_list) = partial_result[previously_assigned_vertex]
                extended_assigned_children = assigned_children + (vertex,)
                for verified_candidate_tuple in verified_candidates_list:
                    (verified_candidate, verified_subsets) = verified_candidate_tuple
                    if verified_candidate not in vertex_ancestors:
                        continue
                    is_correct_assignment = True
                    verified_candidate_access = self.pedigree.vertex_to_ancestor_map[vertex][verified_candidate]
                    new_subsets = dict(verified_subsets)
                    for (image, preimage_size) in verified_subsets.items():
                        image: frozenset
                        new_image = cache_set(image.union(verified_candidate_access))
                        if len(new_image) < preimage_size + 1:
                            is_correct_assignment = False
                            break
                        if preimage_size + 1 > new_subsets.get(new_image, 0):
                            new_subsets[new_image] = preimage_size + 1
                    if not is_correct_assignment:
                        continue
                    candidate_access = cache_set(
                        frozenset(self.pedigree.vertex_to_ancestor_map[vertex][verified_candidate])
                    )
                    if candidate_access not in new_subsets:
                        new_subsets[candidate_access] = 1
                    new_verified_candidates_list.append((verified_candidate, new_subsets))
                if new_verified_candidates_list:
                    verified_assignments = (extended_assigned_children, new_verified_candidates_list)
            if verified_assignments:
                # If there are valid assignments, we can either continue exploring or save the results
                # if this is the last tree vertex
                current_partial_assignment[vertices_coalescent_ids[current_index]] = vertex
                partial_result[vertex] = verified_assignments
                if current_index + 1 < vertices_length:
                    # There are unexplored tree vertices
                    current_index += 1
                    current_stack.append(-1)
                    current_stack.extend(coalescent_vertex_to_candidates[vertices_coalescent_ids[current_index]])
                else:
                    # We have explored all the vertices and verified all the constraints. Now, we can
                    # save the results and switch to the next unexplored branch
                    (assigned_children, candidates_list) = verified_assignments
                    valid_candidates = [(assigned_children, [candidate[0] for candidate in candidates_list])]
                    result.extend(valid_candidates)
                    free_current_vertex_data()
        formatted_result = defaultdict(list)
        for (assigned_children, candidates_list) in result:
            for candidate in candidates_list:
                formatted_result[candidate].append(assigned_children)
        return formatted_result
