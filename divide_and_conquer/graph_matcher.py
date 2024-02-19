"""!
@file graph_matcher.py
@brief This file contains the GraphMatcher class which performs the divide-and-conquer alignment algorithm.
"""

import datetime
import os
import time

from divide_and_conquer.subtree_matcher import *
from genealogical_graph import CoalescentTree
from divide_and_conquer.potential_mrca_processed_graph import *

logs_enabled = True
logs_default_directory_name = "logs"

print_enabled = False


class MatcherLogger:

    def __init__(self, logs_directory_name=logs_default_directory_name):
        if logs_enabled:
            time_in_milliseconds = datetime.datetime.now().timestamp() * 1000
            log_filename = f"{logs_directory_name}/log_{int(time_in_milliseconds)}.txt"
            if not os.path.exists(logs_directory_name):
                os.makedirs(logs_directory_name)
            if print_enabled:
                print(f"The filename is {log_filename}")
            self.file = open(log_filename, 'w')

    def log(self, line: str):
        if logs_enabled:
            self.file.write(f"{line}\n")
            self.file.flush()

    def log_children_domain_space(self, pedigree: PotentialMrcaProcessedGraph,
                                  coalescent_tree: CoalescentTree,
                                  coalescent_tree_vertex: int, coalescent_tree_children: [int],
                                  child_candidate_subtree_matchers_matrix: dict):
        coalescent_tree_level = coalescent_tree.vertex_to_level_map[coalescent_tree_vertex]
        pedigree_level = pedigree.vertex_to_level_map[coalescent_tree_vertex]
        self.log("####################")
        self.log(f"Processing vertex {coalescent_tree_vertex} in the coalescent tree")
        self.log(f"Pedigree level: {pedigree_level}")
        self.log(f"Coalescent tree level: {coalescent_tree_level}")
        self.log(f"There are {len(coalescent_tree.levels) - coalescent_tree_level} vertices above")
        self.log(f"There are {len(pedigree.levels)} levels in the pedigree")
        self.log(f"{coalescent_tree_vertex} has {len(coalescent_tree_children)} children: {coalescent_tree_children}")
        for child in coalescent_tree_children:
            child_candidates = list(child_candidate_subtree_matchers_matrix[child])
            self.log(f"{child} ({len(child_candidate_subtree_matchers_matrix[child])})")
            level_to_candidate = {level: [x for x in child_candidates if pedigree.vertex_to_level_map[x] == level]
                                  for level in range(len(pedigree.levels))}
            for level, vertices in level_to_candidate.items():
                if vertices:
                    self.log(f"{level}: {len(vertices)}")

        self.log("--------------------")

    def log_pair_common_ancestors(self, pedigree: PotentialMrcaProcessedGraph, first_pedigree_vertex: int,
                                  second_pedigree_vertex: int, potential_common_ancestors: [int]):
        self.log(f"{first_pedigree_vertex} ({pedigree.vertex_to_level_map[first_pedigree_vertex]}), "
                 f"{second_pedigree_vertex} ({pedigree.vertex_to_level_map[second_pedigree_vertex]}): "
                 f"{potential_common_ancestors}")

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
        total_candidates_found = 0
        correct_children_assignment_combinations = 0
        for correct_partial_result in inference_result:
            correct_children_assignment_combinations += 1
            (_, valid_candidates) = correct_partial_result
            total_candidates_found += len(valid_candidates)
        self.log(f"Total candidates found: {total_candidates_found}")
        self.log(f"Correct children combinations found: {correct_children_assignment_combinations}")

    def close(self):
        if logs_enabled:
            self.file.close()


class GraphMather:
    """!
    This class runs the divide-and-conquer alignment algorithm on the given preprocessed pedigree and the coalescent
    tree.

    Attributes:
        processed_graph (PotentialMrcaProcessedGraph): The preprocessed pedigree graph that stores the "access vertices"
        through which a descendant vertex can reach the ancestor vertex.
        coalescent_tree (CoalescentTree): The coalescent tree with which the processed_graph is to be aligned.
        logger (MatcherLogger): The logger instance to log steps of the algorithm
        initial_mapping (dict): The initial mapping between the proband vertices in the processed pedigree to
        the vertices in the coalescent tree.
    """

    def __init__(self, processed_graph: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree,
                 logger: MatcherLogger, initial_mapping: dict = None):
        self.pedigree = processed_graph
        self.coalescent_tree = coalescent_tree
        if initial_mapping is None:
            initial_mapping = {x: x for x in processed_graph.probands}
        self.initial_mapping = initial_mapping
        self.logger = logger

    def find_mapping(self):
        """!
        @brief   This method finds all the valid mapping between the given processed pedigree and
                 every sub-clade within the coalescent tree
        @return  Dictionary that maps every vertex in the coalescent tree to the list of valid alignments for
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
        # root_assignments = [coalescent_tree_vertex_to_subtrees[x] for x in self.coalescent_tree.levels[-1]]
        clades = sorted(self.coalescent_tree.get_connected_components(), reverse=True, key=len)[:10]
        if print_enabled:
            print("Filtering")
        start_time = time.time()
        result_dict = dict()
        for clade in clades:
            root_vertex = self.coalescent_tree.get_root_for_clade(clade)
            lowest_non_unary_node = root_vertex
            while len(self.coalescent_tree.children_map[lowest_non_unary_node]) == 1:
                lowest_non_unary_node = self.coalescent_tree.children_map[lowest_non_unary_node][0]
                if lowest_non_unary_node not in self.coalescent_tree.children_map:
                    break
            if lowest_non_unary_node not in self.coalescent_tree.children_map:
                continue
            clade_alignments = []
            for root_matcher in coalescent_tree_vertex_to_subtrees[root_vertex].values():
                root_matcher: SubtreeMatcher
                alignments = root_matcher.get_all_subtree_alignments()
                if print_enabled:
                    print(f"There are {len(alignments)} alignments for the matcher, filtering the results")
                validation_start = time.time()
                valid_alignments = [x for x in alignments if
                                    self.verify_valid_alignment(alignment=x, root_vertex=lowest_non_unary_node)]
                clade_alignments.extend(valid_alignments)
                validation_end = time.time()
                time_taken = validation_end - validation_start
                if print_enabled:
                    print(f"There are {len(valid_alignments)} valid alignments")
                    print(f"Time spent on validation: {time_taken}. Average time taken: {time_taken/len(alignments)}")
            result_dict[root_vertex] = clade_alignments
        end_time = time.time()
        if print_enabled:
            print(f"Time spent on building and filtering the results: {end_time - start_time}")
        return result_dict

    def verify_valid_alignment(self, alignment: dict, root_vertex=None):
        """!
        @brief Verifies that the given alignment is correct
        """
        if root_vertex is None:
            root_vertex = self.coalescent_tree.get_root_for_clade(alignment.keys())
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        target_vertex_label = "t"
        # Adding the edge from the root to the sink vertex
        children_number = len(self.coalescent_tree.children_map[root_vertex])
        assert children_number > 0
        root_vertex_pedigree = alignment[root_vertex]
        network_graph.add_edge(root_vertex_pedigree, target_vertex_label,
                               capacity=children_number)
        proband_number = 0
        for parent in alignment:
            parent_pedigree = alignment[parent]
            if parent not in self.coalescent_tree.children_map:
                network_graph.add_edge(source_vertex_label, parent_pedigree, capacity=1)
                proband_number += 1
                continue
            children = self.coalescent_tree.children_map[parent]
            children_number = len(children)
            assert children_number > 0
            children_pedigree = [alignment[x] for x in children]
            self.pedigree.add_edges_to_mrca_from_descendants(network_graph, parent_pedigree, children_pedigree)
            if parent_pedigree != root_vertex_pedigree:
                network_graph.add_edge(parent_pedigree, target_vertex_label, capacity=children_number - 1)
        maximum_flow = networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=target_vertex_label)
        return maximum_flow == proband_number

    def get_subtrees_from_children(self, focal_vertex: int, vertex_subtree_dict: {int: {int: SubtreeMatcher}}):
        """!
        @brief This method finds all the valid alignments for the sub-clade where the given focal_vertex is the root.
        @param focal_vertex The root of the sub-clade for which the inference is to be performed
        @param vertex_subtree_dict The dictionary containing all the valid alignment for the focal vertex
        @return The list of all the valid sub-clade alignments for the specified focal vertex. The resulting
                list consists of \link SubtreeMatcher subtree_matcher.SubtreeMatcher \endlink objects
        """
        focal_vertex_children: [int] = self.coalescent_tree.children_map[focal_vertex]
        if not focal_vertex_children:
            raise Exception("Isolated vertex in the coalescent tree")
        if len(focal_vertex_children) == 1:
            child, = focal_vertex_children
            result = {x: SubtreeMatcher(root_coalescent_tree=focal_vertex, root_pedigree=x,
                                        subtrees_matchers=[{child: y}]) for x, y in vertex_subtree_dict[child].items()}
            return result
        if print_enabled:
            print("####################")
            print(f"Inference for {focal_vertex} ({self.coalescent_tree.vertex_to_level_map[focal_vertex]}),"
                  f" there are {len(focal_vertex_children)} children")
            for child in focal_vertex_children:
                print(f"{child}: {len(vertex_subtree_dict[child])}")
        inference_start = time.time()
        (children_order, inference_result) = self.pedigree.get_pmracs_for_vertices(
            coalescent_vertex_to_candidates=vertex_subtree_dict,
            vertices_coalescent_ids=focal_vertex_children
        )
        inference_end = time.time()
        time_taken = inference_end - inference_start
        if print_enabled:
            print(f"Inference time: {time_taken}")
            print("####################")
        self.logger.log_vertex_inference_time(self.coalescent_tree, time_taken, focal_vertex, focal_vertex_children,
                                              vertex_subtree_dict, inference_result)
        result_build_start = time.time()
        candidate_subtree_matcher_dictionary = dict()
        for single_result in inference_result:
            (assigned_children, focal_vertex_candidates) = single_result
            assert len(children_order) == len(assigned_children) == len(focal_vertex_children)
            children_dictionary = dict()
            for index, child_assigned_pedigree_id in enumerate(assigned_children):
                child_coalescent_id = children_order[index]
                assigned_child_matcher = vertex_subtree_dict[child_coalescent_id][child_assigned_pedigree_id]
                children_dictionary[child_coalescent_id] = assigned_child_matcher
            for focal_vertex_candidate in focal_vertex_candidates:
                subtree_matcher = candidate_subtree_matcher_dictionary.setdefault(focal_vertex_candidate,
                                                                                  SubtreeMatcher(focal_vertex,
                                                                                                 focal_vertex_candidate,
                                                                                                 []))
                subtree_matcher.children_assignments.append(children_dictionary)
        result_build_end = time.time()
        if print_enabled:
            print(f"Building the result {result_build_end - result_build_start}")
            print(f"There are {len(candidate_subtree_matcher_dictionary)} resulting assignments")
        return candidate_subtree_matcher_dictionary


def get_subtree_matcher_for_coalescent_tree_proband(proband: int, proband_pedigree_id: int):
    """
    @brief Helper function returning the valid sub-alignments for a proband (which is simply one identity alignment).
    @param proband The id of the proband in the coalescent tree
    @param proband_pedigree_id The id of the corresponding vertex in the pedigree
    """
    subtree_matcher = SubtreeMatcher(root_coalescent_tree=proband, root_pedigree=proband_pedigree_id)
    return {proband_pedigree_id: subtree_matcher}
