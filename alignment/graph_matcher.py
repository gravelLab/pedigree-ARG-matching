"""!
@file graph_matcher.py
@brief This file contains the GraphMatcher class which performs the divide-and-conquer alignment algorithm.
"""

import datetime
import os
import time
from collections import defaultdict
from dataclasses import dataclass
from functools import reduce
from pathlib import Path

import networkx

from alignment.configuration import *
from alignment.potential_mrca_processed_graph import *
from alignment.subtree_matcher import *
from graph.coalescent_tree import CoalescentTree
from alignment.log import print_log
from scripts.utility import dict_has_duplicate_values
import yaml


class MatcherLogger:

    def __init__(self, logs_directory_path=logs_default_directory_name):
        time_in_milliseconds = datetime.datetime.now().timestamp() * 1000
        log_filename = Path(logs_directory_path) / f"log_{int(time_in_milliseconds)}.txt"
        if not os.path.exists(logs_directory_path):
            os.makedirs(logs_directory_path)
        print_log(f"The filename is {log_filename}")
        self.file = open(log_filename, 'w')

    def log(self, line: str):
        self.file.write(f"{line}\n")
        self.file.flush()

    def log_vertex_inference_time(self, coalescent_tree: CoalescentTree,
                                  spent_time: float, focal_vertex_coalescent_id: int, focal_vertex_children: [int],
                                  child_candidates_map: {int: [int]}, inference_result):
        self.log(f"Vertex level: {coalescent_tree.vertex_to_level_map[focal_vertex_coalescent_id]}")
        self.log(f"Solving the inference for {focal_vertex_coalescent_id} which has {len(focal_vertex_children)} "
                 f"children")
        self.log(f"The children's domain:")
        for child in focal_vertex_children:
            self.log(f"{child}: {len(child_candidates_map[child])}")
        self.log(f"Time taken for vertex inference: {spent_time}")
        total_candidates_found = len(inference_result)
        self.log(f"Total candidates found: {total_candidates_found}")

    def close(self):
        self.file.close()


def get_subtree_matcher_for_coalescent_tree_proband(proband: int, proband_pedigree_ids: [int]):
    """
    @brief Helper function returning the valid sub-alignments for a proband (which is simply one identity alignment).
    @param proband The id of the proband in the coalescent tree
    @param proband_pedigree_ids The ids of the corresponding vertices in the pedigree
    """
    return {proband_pedigree_id: SubtreeMatcher(root_coalescent_tree=proband, root_pedigree=proband_pedigree_id)
            for proband_pedigree_id in proband_pedigree_ids}


def get_initial_mapping_for_mode(coalescent_tree: CoalescentTree, mode: InitialMatchingMode = current_matching_mode):
    if mode == InitialMatchingMode.PLOID:
        return {x: [x] for x in coalescent_tree.probands}
    return {x: [2 * (x // 2), (2 * (x // 2) + 1)] for x in coalescent_tree.probands}


class YAMLValidationError(Exception):
    """Custom exception for YAML validation errors."""
    pass


class PloidType(Enum):
    Paternal = "P"
    Maternal = "M"
    Unphased = "U"


@dataclass
class InitialAssignment:
    coalescent_id: int
    pedigree_id: int
    ploid_type: PloidType


def validate_and_parse_yaml(file_path: str) -> dict[int, InitialAssignment]:
    """
    Validates the specified YAML file and parses it into a dictionary mapping
    coalescent_id to an InitialAssignment object.

    :param file_path: Path to the YAML file.
    :return: Dictionary where keys are coalescent_ids and values are InitialAssignment objects.
    :raises: YAMLValidationError if the file is invalid.
    """
    required_keys = {"initial_assignments", "coalescent_tree_path", "pedigree_path"}

    try:
        # Load the YAML file
        with open(file_path, 'r') as f:
            data = yaml.safe_load(f)

        # Check if all required keys are present
        if not required_keys.issubset(data.keys()):
            raise YAMLValidationError(f"Missing required keys. Expected: {required_keys}, Found: {data.keys()}")

        # Validate "initial_assignments" format
        if not isinstance(data["initial_assignments"], list):
            raise YAMLValidationError("initial_assignments should be a list.")

        coalescent_dict = {}

        for entry in data["initial_assignments"]:
            if not isinstance(entry, dict):
                raise YAMLValidationError(f"Each entry in initial_assignments must be a dictionary. Found: {entry}")
            if not {"coalescent_id", "pedigree_id", "ploid_type"}.issubset(entry.keys()):
                raise YAMLValidationError(f"Missing keys in initial_assignments entry. Found: {entry.keys()}")

            # Validate and convert ploid_type
            try:
                ploid_type = PloidType(entry["ploid_type"])
            except ValueError:
                raise YAMLValidationError(
                    f"Invalid ploid_type: {entry['ploid_type']}. Must be one of {[e.value for e in PloidType]}"
                )

            # Validate coalescent_id uniqueness
            coalescent_id = int(entry["coalescent_id"])

            # Check if coalescent_id already exists
            if coalescent_id in coalescent_dict:
                raise YAMLValidationError(f"Duplicate coalescent_id found: {coalescent_id}")

            # Create InitialAssignment object
            initial_assignment = InitialAssignment(
                coalescent_id=coalescent_id,
                pedigree_id=int(entry["pedigree_id"]),
                ploid_type=ploid_type
            )

            # Add to the dictionary
            coalescent_dict[coalescent_id] = initial_assignment

        # Validate paths
        if not os.path.isfile(data["coalescent_tree_path"]):
            raise YAMLValidationError(f"Invalid coalescent_tree_path: {data['coalescent_tree_path']}")
        if not os.path.isfile(data["pedigree_path"]):
            raise YAMLValidationError(f"Invalid pedigree_path: {data['pedigree_path']}")

        return coalescent_dict

    except FileNotFoundError:
        raise YAMLValidationError(f"File not found: {file_path}")
    except yaml.YAMLError as e:
        raise YAMLValidationError(f"Error parsing YAML file: {e}")
    except Exception as e:
        raise YAMLValidationError(f"Unexpected error: {e}")


class GraphMatcher:
    """!
    This class runs the divide-and-conquer alignment algorithm on the given preprocessed pedigree and the coalescent
    tree.
    """

    def __init__(self, processed_graph: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree,
                 matching_mode: MatchingMode = None,
                 logger: MatcherLogger = None, initial_mapping: dict = None):
        """
        @param processed_graph (PotentialMrcaProcessedGraph): The preprocessed pedigree graph that stores the
        "access vertices" through which a descendant vertex can reach the ancestor vertex.
        @param coalescent_tree (CoalescentTree): The coalescent tree with which the processed_graph is to be aligned.
        @param logger (MatcherLogger): The logger instance to log steps of the algorithm
        @param initial_mapping (dict): The initial mapping between the proband vertices in the processed pedigree to
        the vertices in the coalescent tree.
        """
        self.pedigree = processed_graph
        self.coalescent_tree = coalescent_tree
        if matching_mode is None:
            matching_mode = current_matching_mode
        self.matching_mode = matching_mode
        if initial_mapping is None:
            initial_mapping = get_initial_mapping_for_mode(coalescent_tree=self.coalescent_tree)
        self.initial_mapping = initial_mapping
        self.logger = logger

    def log(self, string: str):
        if self.logger:
            self.logger.log(string)

    def find_mapping(self):
        """!
        @brief   This method finds all the valid mapping between the given processed pedigree and
                 every sub-clade within the coalescent tree
        @return  Dictionary that maps every clade root vertex in the coalescent tree to the list of valid alignments for
                 the sub-clade where the chosen vertex is the root.
                 The alignments are represented as a list of SubtreeMatcher objects
        """
        coalescent_tree_vertex_to_subtrees = dict()
        for vertex in self.coalescent_tree.levels[0]:
            coalescent_tree_vertex_to_subtrees[vertex] = get_subtree_matcher_for_coalescent_tree_proband(
                vertex, self.initial_mapping[vertex])
        for level in self.coalescent_tree.levels[1:]:
            for vertex in level:
                coalescent_tree_vertex_to_subtrees[vertex] = self.get_subtrees_from_children(
                    vertex,
                    coalescent_tree_vertex_to_subtrees)
        # Filtering the map, so that only the assignments for the roots of the clades remain
        print_log("Filtering the alignments")
        coalescent_tree_vertex_to_subtrees = {
            x: value
            for x, value in coalescent_tree_vertex_to_subtrees.items()
            if x not in self.coalescent_tree.parents_map
        }
        start_time = time.time()
        result_dict = self.filter_alignments(coalescent_tree_vertex_to_subtrees)
        end_time = time.time()
        print_log(f"Time spent on building and filtering the results: {end_time - start_time}")
        return result_dict

    def filter_alignments(self, coalescent_tree_vertex_to_subtrees: dict):
        if self.matching_mode == MatchingMode.ALL_ALIGNMENTS:
            # return self.filter_alignments_iteratively_recursive_cache(coalescent_tree_vertex_to_subtrees)
            return self.filter_alignments_brute_force(coalescent_tree_vertex_to_subtrees)
        return self.filter_alignments_example_per_root_assignment(coalescent_tree_vertex_to_subtrees)

    def filter_alignments_example_per_root_assignment(self, coalescent_tree_vertex_to_subtrees: dict):
        result_dict = dict()
        for clade_root_vertex, clade_root_vertex_assignments in coalescent_tree_vertex_to_subtrees.items():
            clade_alignments = []
            root_vertex_matchers: [SubtreeMatcher] = clade_root_vertex_assignments.values()
            self.log(f"Looking for valid assignments for the root vertex {clade_root_vertex}"
                     f" There are {len(root_vertex_matchers)}"
                     f"unverified assignments to the root of the clade")
            for root_matcher in root_vertex_matchers:
                assert clade_root_vertex == root_matcher.root_coalescent_tree
                self.log(f"Finding an example alignment for the assignment "
                         f"{clade_root_vertex}: {root_matcher.root_pedigree}")
                alignments = root_matcher.get_all_subtree_alignments()
                self.log(f"There are {len(alignments)} unverified alignments")
                for alignment in alignments:
                    if self.verify_valid_alignment(alignment=alignment,
                                                   root_vertex=root_matcher.root_coalescent_tree):
                        clade_alignments.append(alignment)
                        break
            result_dict[clade_root_vertex] = clade_alignments
        return result_dict

    def filter_clade_root_assignments_independently(self, coalescent_tree_vertex_to_subtrees: dict,
                                                    matcher_filter_function):
        result_dict = dict()
        for clade_root_vertex, clade_root_vertex_assignments in coalescent_tree_vertex_to_subtrees.items():
            clade_alignments = []
            root_vertex_matchers: [SubtreeMatcher] = clade_root_vertex_assignments.values()
            self.log(f"Filtering the alignments, there are {len(root_vertex_matchers)}"
                     f"unverified assignments to the root of the clade")
            for root_matcher in root_vertex_matchers:
                self.log(f"Filtering the alignments for the root matcher "
                         f"{root_matcher.root_pedigree}")
                validation_start = time.time()
                valid_alignments = matcher_filter_function(root_matcher)
                validation_end = time.time()
                clade_alignments.extend(valid_alignments)
                time_taken = validation_end - validation_start
                print_log(f"There are {len(valid_alignments)} valid alignments")
                print_log(f"Time spent on validation: {time_taken}.")
                clade_alignments.extend(valid_alignments)
            result_dict[clade_root_vertex] = clade_alignments
        return result_dict

    def filter_alignments_brute_force(self, coalescent_tree_vertex_to_subtrees: dict):
        return self.filter_clade_root_assignments_independently(
            coalescent_tree_vertex_to_subtrees=coalescent_tree_vertex_to_subtrees,
            matcher_filter_function=self.get_verified_alignments_brute_force)

    def filter_alignments_product_recursive(self, coalescent_tree_vertex_to_subtrees: dict):
        return self.filter_clade_root_assignments_independently(
            coalescent_tree_vertex_to_subtrees=coalescent_tree_vertex_to_subtrees,
            matcher_filter_function=self.get_verified_alignments_product_recursive)

    def filter_alignments_iteratively_recursive_cache(self, coalescent_tree_vertex_to_subtrees: dict):
        # TODO: Make a separate non-cache implementation
        result_dict = dict()
        sub_alignments_cache = dict()
        for clade_root_vertex, clade_root_vertex_assignments in coalescent_tree_vertex_to_subtrees.items():
            clade_alignments = []
            root_vertex_matchers: [SubtreeMatcher] = clade_root_vertex_assignments.values()
            self.log(f"Filtering the alignments, there are {len(root_vertex_matchers)}"
                     f"unverified assignments to the root of the clade")
            for root_matcher in root_vertex_matchers:
                self.log(f"Filtering the alignments for the root matcher "
                         f"{root_matcher.root_pedigree}")
                validation_start = time.time()
                valid_alignments = self.get_verified_alignments_iteratively_recursive_cache(
                    matcher=root_matcher,
                    sub_alignments_cache=sub_alignments_cache)
                validation_end = time.time()
                clade_alignments.extend(valid_alignments)
                time_taken = validation_end - validation_start
                print_log(f"There are {len(valid_alignments)} valid alignments")
                print_log(f"Time spent on validation: {time_taken}.")
                clade_alignments.extend(valid_alignments)
            result_dict[clade_root_vertex] = clade_alignments
        return result_dict

    def get_verified_alignments_product_recursive(self, matcher: SubtreeMatcher):
        if matcher.subtree_alignments is not None:
            return matcher.subtree_alignments
        results = []
        root_assignment_dict = {matcher.root_coalescent_tree: matcher.root_pedigree}
        # If the given vertex has no children, the resulting alignment is a dictionary with one key-value pair
        if matcher.children_assignments is None:
            return [root_assignment_dict]
        # If there are children assignments, loop over all the corresponding children alignments
        for children_assignment in matcher.children_assignments:
            children_alignments = [self.get_verified_alignments_product_recursive(matcher=x) for x in
                                   children_assignment.values()]
            search_space = reduce(lambda x, y: x * len(y), children_alignments, 1)
            print_log(f"Search space {search_space}")
            for children_dictionaries in product(*children_alignments):
                new_result = dict(root_assignment_dict)
                for dictionary in children_dictionaries:
                    new_result.update(dictionary)
                if self.verify_valid_alignment(alignment=new_result, root_vertex=matcher.root_coalescent_tree):
                    results.append(new_result)
        matcher.subtree_alignments = results
        return results

    def get_verified_alignments_iteratively_recursive_cache(self, matcher: SubtreeMatcher,
                                                            sub_alignments_cache: dict,
                                                            print_debug=True):
        #  When the number of unverified alignments is quite small, this approach can be slower
        if matcher.subtree_alignments is not None:
            return matcher.subtree_alignments
        results = []
        root_assignment_dict = {matcher.root_coalescent_tree: matcher.root_pedigree}
        # If the given vertex has no children, the resulting alignment is a dictionary with one key-value pair
        if matcher.children_assignments is None:
            return [root_assignment_dict]
        # There are children assignments, loop over all the corresponding children alignments
        for index, children_assignment in enumerate(matcher.children_assignments):
            if print_debug:
                print(f"Processing assignment {index + 1} of {len(matcher.children_assignments)}")
            # current_valid_alignments = [root_assignment_dict]
            children_matchers = sorted(children_assignment.values(),
                                       key=lambda x: x.root_coalescent_tree in self.coalescent_tree.children_map,
                                       reverse=False)
            current_valid_alignments = self.get_verified_alignments_iteratively_recursive_cache(
                children_matchers[0], sub_alignments_cache)
            for child_index, child_matcher in enumerate(children_matchers[1:]):
                if print_debug:
                    print(f"Processing child {child_index + 2} "
                          f"of {len(self.coalescent_tree.children_map[matcher.root_coalescent_tree])}")
                    print(f"Current number of sub-alignments is {len(current_valid_alignments)}")
                current_children_matchers = children_matchers[:child_index + 2]
                children_assignment = tuple(
                    (child_matcher.root_coalescent_tree, child_matcher.root_pedigree)
                    for child_matcher in current_children_matchers
                )
                new_valid_alignments = sub_alignments_cache.get(children_assignment, [])
                # new_valid_alignments = []
                if new_valid_alignments:
                    print("Using cached results")
                else:
                    print("#########################")
                    child_matcher_alignments = self.get_verified_alignments_iteratively_recursive_cache(
                        matcher=child_matcher,
                        sub_alignments_cache=sub_alignments_cache)
                    print("#########################")
                    for child_matcher_alignment in child_matcher_alignments:
                        for current_valid_alignment in current_valid_alignments:
                            new_alignment = dict(child_matcher_alignment)
                            new_alignment.update(current_valid_alignment)
                            if self.verify_valid_alignment_for_potential_root(alignment=new_alignment):
                                new_valid_alignments.append(new_alignment)
                            else:
                                print("Discard")
                    sub_alignments_cache[children_assignment] = new_valid_alignments
                current_valid_alignments = new_valid_alignments
            result_alignments = []
            for current_valid_alignment in current_valid_alignments:
                new_alignment = dict(current_valid_alignment)
                new_alignment.update(root_assignment_dict)
                if self.verify_valid_alignment_for_potential_root(alignment=new_alignment):
                    result_alignments.append(new_alignment)
            results.extend(result_alignments)
        matcher.subtree_alignments = results
        return results

    def get_verified_alignments_brute_force(self, matcher: SubtreeMatcher):
        alignments = matcher.get_all_subtree_alignments()
        print_log(f"There are {len(alignments)} alignments for the matcher, filtering the results")
        # valid_alignments = [x for x in alignments if
        #                     self.verify_valid_alignment(alignment=x, root_vertex=matcher.root_coalescent_tree)]
        valid_alignments = [x for x in alignments if
                            self.verify_valid_alignment_for_potential_root(alignment=x)]
        return valid_alignments

    def verify_valid_alignment_for_potential_root(self, alignment: dict):
        if dict_has_duplicate_values(alignment):
            return False
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        target_vertex_label = "t"
        proband_number = 0
        for tree_vertex in alignment:
            tree_vertex_pedigree = alignment[tree_vertex]
            tree_vertex_parents = [x for x in self.coalescent_tree.parents_map.get(tree_vertex, []) if x in alignment]
            if tree_vertex not in self.coalescent_tree.children_map:
                network_graph.add_edge(source_vertex_label, tree_vertex_pedigree, capacity=1)
                proband_number += 1
                if not tree_vertex_parents:
                    network_graph.add_edge(tree_vertex_pedigree, target_vertex_label, capacity=1)
                continue
            children = self.coalescent_tree.children_map[tree_vertex]
            assert len(children) > 0
            children_pedigree = [alignment[x] for x in children if alignment.get(x) is not None]
            # The number of children present in the alignment
            children_number = len(children_pedigree)
            self.add_edges_to_mrca_from_descendants(network_graph, tree_vertex_pedigree, children_pedigree)
            if not tree_vertex_parents:
                network_graph.add_edge(tree_vertex_pedigree, target_vertex_label, capacity=children_number)
            else:
                network_graph.add_edge(tree_vertex_pedigree, target_vertex_label, capacity=children_number - 1)
        maximum_flow = networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=target_vertex_label)
        return maximum_flow == proband_number

    def verify_valid_alignment(self, alignment: dict, root_vertex: int):
        """!
        @brief Verifies that the given alignment is correct
        """
        if dict_has_duplicate_values(alignment):
            return False
        # if root_vertex is None:
        #     root_vertex = self.coalescent_tree.get_root_for_clade(alignment.keys())
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        target_vertex_label = "t"
        # Adding the edge from the root to the sink vertex
        root_vertex_children_number = len(self.coalescent_tree.children_map[root_vertex])
        assert root_vertex_children_number > 0
        root_vertex_pedigree = alignment[root_vertex]
        network_graph.add_edge(root_vertex_pedigree, target_vertex_label,
                               capacity=root_vertex_children_number)
        proband_number = 0
        for parent in alignment:
            parent_pedigree = alignment[parent]
            if parent not in self.coalescent_tree.children_map:
                network_graph.add_edge(source_vertex_label, parent_pedigree, capacity=1)
                proband_number += 1
                continue
            children = self.coalescent_tree.children_map[parent]
            assert len(children) > 0
            children_pedigree = [alignment[x] for x in children if alignment.get(x) is not None]
            # if len(children_pedigree) != len(children):
            #     print("WARNING: not all the children are present in the mapping")
            children_number = len(children_pedigree)
            self.add_edges_to_mrca_from_descendants(network_graph, parent_pedigree, children_pedigree)
            if parent_pedigree != root_vertex_pedigree:
                network_graph.add_edge(parent_pedigree, target_vertex_label, capacity=children_number - 1)
        maximum_flow = networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=target_vertex_label)
        return maximum_flow == proband_number

    def add_edges_to_mrca_from_descendants(self, flow_network: networkx.DiGraph, mrca: int, descendant_vertices: [int]):
        """
        :param flow_network:
        :param mrca:
        :param descendant_vertices:
        :return:
        """
        mrca_access_label = f"{mrca}'"
        flow_network.add_edge(mrca_access_label, mrca, capacity=len(descendant_vertices))
        for descendant in descendant_vertices:
            # descendant_mrca_access_vertices = self.vertex_to_ancestor_map[descendant][mrca]
            # flow_network.add_edges_from([(x, mrca_access_label) for x in descendant_mrca_access_vertices], capacity=1)
            # current_level_vertices = descendant_mrca_access_vertices
            current_level_vertices = [mrca]
            while current_level_vertices:
                next_level_vertices = set()
                for vertex in current_level_vertices:
                    if vertex == descendant:
                        continue
                    neighbors = self.pedigree.vertex_to_ancestor_map[descendant][vertex]
                    for neighbor in neighbors:
                        vertex_access_label = f"{vertex}'"
                        # We can potentially override the edge's capacity if vertex is the mrca in the
                        # coalescent tree
                        if not flow_network.has_edge(vertex_access_label, vertex):
                            flow_network.add_edge(vertex_access_label, vertex, capacity=1)
                        flow_network.add_edge(neighbor, vertex_access_label, capacity=1)
                    next_level_vertices.update(neighbors)
                current_level_vertices = next_level_vertices

    def get_subtrees_from_children(self, focal_vertex: int, vertex_subtree_dict: {int: {int: SubtreeMatcher}}):
        """!
        @brief This method finds all the valid alignments for the sub-clade where the given focal_vertex is the root.
        @param focal_vertex The root of the sub-clade for which the inference is to be performed
        @param vertex_subtree_dict The dictionary containing all the valid alignment for the focal vertex
        @return The list of all the valid sub-clade alignments for the specified focal vertex. The resulting
                list consists of SubtreeMatcher subtree_matcher.SubtreeMatcher objects
        """
        focal_vertex_children: [int] = self.coalescent_tree.children_map[focal_vertex]
        if not focal_vertex_children:
            raise Exception("Isolated vertex in the coalescent tree")
        if len(focal_vertex_children) == 1:
            child, = focal_vertex_children
            result = {x: SubtreeMatcher(root_coalescent_tree=focal_vertex, root_pedigree=x,
                                        children_assignments=[{child: y}]) for x, y in
                      vertex_subtree_dict[child].items()}
            return result
        if print_enabled:
            print("####################")
            print(f"Inference for {focal_vertex} ({self.coalescent_tree.vertex_to_level_map[focal_vertex]}),"
                  f" there are {len(focal_vertex_children)} children")
            for child in focal_vertex_children:
                print(f"{child}: {len(vertex_subtree_dict[child])}")
        inference_start = time.time()
        (children_order, inference_result) = self.get_pmracs_for_vertices(
            coalescent_vertex_to_candidates=vertex_subtree_dict,
            vertices_coalescent_ids=focal_vertex_children
        )
        inference_end = time.time()
        time_taken = inference_end - inference_start
        print_log(f"Inference time: {time_taken}")
        print_log("####################")
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
        # for single_result in inference_result:
        #     (assigned_children, focal_vertex_candidates) = single_result
        #     assert len(children_order) == len(assigned_children) == len(focal_vertex_children)
        #     children_dictionary = dict()
        #     for index, child_assigned_pedigree_id in enumerate(assigned_children):
        #         child_coalescent_id = children_order[index]
        #         assigned_child_matcher = vertex_subtree_dict[child_coalescent_id][child_assigned_pedigree_id]
        #         children_dictionary[child_coalescent_id] = assigned_child_matcher
        #     for focal_vertex_candidate in focal_vertex_candidates:
        #         subtree_matcher = candidate_subtree_matcher_dictionary.setdefault(focal_vertex_candidate,
        #                                                                           SubtreeMatcher(focal_vertex,
        #                                                                                          focal_vertex_candidate,
        #                                                                                          []))
        #         subtree_matcher.children_assignments.append(children_dictionary)
        result_build_end = time.time()
        print_log(f"Building the result {result_build_end - result_build_start}")
        print_log(f"There are {len(candidate_subtree_matcher_dictionary)} resulting assignments")
        return candidate_subtree_matcher_dictionary

    # ----------------------------------------- Alignment logic ------------------------------------------------------

    def verify_pmrca_for_vertex_pair(self, first_vertex: int, second_vertex: int, common_ancestor: int):
        """!
        @brief This function verifies that the given vertex is a potential mrca for the given pair of vertices.
        @param first_vertex The first vertex from the pair.
        @param second_vertex The second vertex from the pair.
        @param common_ancestor A common ancestor of the given pair of vertices that is being verified.
        @return True if the given vertex is potential mrca for the given pair of vertices, False otherwise.
        """
        first_vertex_access_vertices = self.pedigree.vertex_to_ancestor_map[first_vertex][common_ancestor]
        return not (len(first_vertex_access_vertices) == 1 and
                    first_vertex_access_vertices ==
                    self.pedigree.vertex_to_ancestor_map[second_vertex][common_ancestor])

    def triplet_condition_holds(self, first_vertex: int, second_vertex: int, third_vertex: int,
                                potential_mrca_candidate: int):
        """!
        @brief This function verifies that the potential mrca candidate satisfies all the conditions involving the
        third vertex. More specifically, it verifies that the
        <a href="https://en.wikipedia.org/wiki/Hall%27s_marriage_theorem">Hall's condition</a> is satisfied for the
        triplet (first_vertex, second_vertex, third_vertex).
        @param first_vertex: The first vertex from the triplet.
        @param second_vertex: The second vertex from the triplet.
        @param third_vertex: The third vertex from the triplet.
        @param potential_mrca_candidate: A common ancestor of the given triplet of vertices that satisfies the
        Hall's condition for the (first_vertex, second_vertex) pair.
        @return True if potential mrca candidate is indeed a potential mrca, False otherwise.
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

    def filter_common_ancestors_for_vertex_pair(self, first_candidate: int, second_candidate: int):
        """!
        @brief Finds the pmracs vertices for the given vertex pair. Important: This method does not
        check whether the passed vertices are different
        @param first_candidate: The first vertex id
        @param second_candidate: The second vertex id
        @return: The list of pmracs for the given vertices
        """
        verified_ancestors = []
        if (not self.pedigree.parents_map.get(first_candidate, []) or
                not self.pedigree.parents_map.get(second_candidate, [])):
            return verified_ancestors
        while self.pedigree.parents_map[first_candidate] == self.pedigree.parents_map[second_candidate]:
            verified_ancestors.extend(self.pedigree.parents_map[first_candidate])
            [first_candidate, second_candidate] = self.pedigree.parents_map[first_candidate]
            if (not self.pedigree.parents_map.get(first_candidate, []) or
                    not self.pedigree.parents_map.get(second_candidate, [])):
                return verified_ancestors
        first_ancestors = self.pedigree.vertex_to_ancestor_map[first_candidate]
        second_ancestors = self.pedigree.vertex_to_ancestor_map[second_candidate]
        if len(second_ancestors) < len(first_ancestors):
            first_ancestors, second_ancestors = second_ancestors, first_ancestors
        verified_ancestors.extend([ancestor for ancestor in first_ancestors if ancestor in second_ancestors and
                                   self.verify_pmrca_for_vertex_pair(first_candidate, second_candidate, ancestor)
                                   ])
        return verified_ancestors

    def get_pmracs_for_vertex_pair(self, first: int, second: int, coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief Explores the possible candidates for the given coalescent vertices and finds the corresponding pmracs
        @param first: The first coalescent vertex id
        @param second: The second coalescent vertex id
        @param coalescent_vertex_to_candidates: A dictionary mapping coalescent vertex ids to the list of pedigree
        vertices
        @return The dictionary mapping a valid pmrca id to the corresponding pedigree assignments for the passed
        coalescent vertices
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
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(
                    first_vertex_candidate,
                    second_vertex_candidate)
                for verified_ancestor in verified_ancestors:
                    result[verified_ancestor].append((first_vertex_candidate, second_vertex_candidate))
        return result

    def get_pmracs_for_vertex_triple_iterative(self, first: int, second: int, third: int,
                                               coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief Explores the possible candidates for the given coalescent vertices and finds the corresponding pmracs
        @param first: The first coalescent vertex
        @param second:The second coalescent vertex
        @param third: The third coalescent vertex
        @param coalescent_vertex_to_candidates: A dictionary mapping coalescent vertex ids to the list of pedigree
        vertices
        @return The dictionary mapping a valid pmrca id to the corresponding pedigree assignments for the passed
        coalescent vertices
        """
        triplet_cache = dict()

        def get_triplet_tuple(first_element, second_element, third_element):
            triplet_tuple = (first_element, second_element, third_element)
            if triplet_tuple in triplet_cache:
                return triplet_cache[triplet_tuple]
            triplet_cache[triplet_tuple] = triplet_tuple
            return triplet_tuple

        first_second_pair_result = self.get_pmracs_for_vertex_pair(first, second, coalescent_vertex_to_candidates)
        first_third_pair_result = self.get_pmracs_for_vertex_pair(first, third, coalescent_vertex_to_candidates)
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
                        if (self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, shared_common_ancestor)
                                and self.triplet_condition_holds(first_candidate, second_candidate,
                                                                 third_candidate, shared_common_ancestor)):
                            result[shared_common_ancestor].append(current_triplet_tuple)
        return result

    def get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
                                coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief This function calculates the potential most recent common ancestor (pmrca) for the given vertices.
        @param vertices_coalescent_ids The ids for the coalescent vertices.
        @param coalescent_vertex_to_candidates: A dictionary mapping an id of a coalescent vertex to the list of
        ids of pedigree vertices that can be assigned to this coalescent vertex.
        @return All the valid pmracs for the given vertices.
        """
        vertices_length = len(vertices_coalescent_ids)
        vertices_coalescent_ids = sorted(vertices_coalescent_ids,
                                         key=lambda child: len(coalescent_vertex_to_candidates[child]),
                                         reverse=False
                                         )
        if vertices_length == 2:
            result = self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                     coalescent_vertex_to_candidates)
        elif vertices_length == 3:
            result = self.get_pmracs_for_vertex_triple_iterative(vertices_coalescent_ids[0],
                                                                 vertices_coalescent_ids[1],
                                                                 vertices_coalescent_ids[2],
                                                                 coalescent_vertex_to_candidates)
        else:
            result = self.get_pmracs_for_vertices_dfs(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        return vertices_coalescent_ids, result

    def get_pmracs_for_vertices_dfs(self, vertices_coalescent_ids: [int],
                                    coalescent_vertex_to_candidates: {int: [int]}):
        """
        Explores the search graph in a dfs manner storing the results only for the currently explored branch
        @param vertices_coalescent_ids: The ids for the coalescent vertices.
        @param coalescent_vertex_to_candidates: A dictionary mapping an id of a coalescent vertex to the list of
        ids of pedigree vertices that can be assigned to this coalescent vertex.
        @return: All the valid pmracs for the given vertices.
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
        current_stack = []
        current_stack.extend(coalescent_vertex_to_candidates[vertices_coalescent_ids[0]])
        set_cache = dict()

        def cache_set(set_to_cache: frozenset):
            if set_to_cache in set_cache:
                return set_cache[set_to_cache]
            set_cache[set_to_cache] = set_to_cache
            return set_to_cache

        def free_current_vertex_data():
            nonlocal current_index
            current_coalescent_vertex = vertices_coalescent_ids[current_index]
            current_fixed_vertex = current_partial_assignment[current_coalescent_vertex]
            partial_result.pop(current_fixed_vertex)
            current_partial_assignment.pop(current_coalescent_vertex)

        while current_stack:
            vertex = current_stack.pop()
            if vertex in partial_result:
                continue
            if vertex == -1:
                current_index -= 1
                free_current_vertex_data()
                continue
            vertex_ancestors = self.pedigree.get_vertex_ancestors(vertex)
            verified_assignments = None
            if not current_index:
                first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(vertex)
                                             if len(self.pedigree.children_map[x]) >= vertices_length]
                if first_candidate_ancestors:
                    verified_assignments = (
                        (vertex,),
                        ([(ancestor, {cache_set(frozenset(self.pedigree.vertex_to_ancestor_map[vertex][ancestor])): 1})
                          for ancestor in vertex_ancestors])
                    )
            else:
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
                current_partial_assignment[vertices_coalescent_ids[current_index]] = vertex
                partial_result[vertex] = verified_assignments
                if current_index + 1 < vertices_length:
                    current_index += 1
                    current_stack.append(-1)
                    current_stack.extend(coalescent_vertex_to_candidates[vertices_coalescent_ids[current_index]])
                else:
                    (assigned_children, candidates_list) = verified_assignments
                    valid_candidates = [(assigned_children, [candidate[0] for candidate in candidates_list])]
                    result.extend(valid_candidates)
                    free_current_vertex_data()
        formatted_result = defaultdict(list)
        for (assigned_children, candidates_list) in result:
            for candidate in candidates_list:
                formatted_result[candidate].append(assigned_children)
        return formatted_result
