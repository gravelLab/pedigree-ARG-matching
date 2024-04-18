"""!
@file graph_matcher.py
@brief This file contains the GraphMatcher class which performs the divide-and-conquer alignment algorithm.
"""

import datetime
import os
from collections import defaultdict
from enum import Enum

from diskcache import Cache

from divide_and_conquer.potential_mrca_processed_graph import *
from divide_and_conquer.subtree_matcher import *
from genealogical_graph import CoalescentTree

logs_enabled = True
logs_default_directory_name = "logs"

print_enabled = True


class InitialMatchingMode(Enum):
    PLOID = 1
    INDIVIDUAL = 2


current_initial_matching_mode = InitialMatchingMode.INDIVIDUAL


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
        total_candidates_found = len(inference_result)
        # correct_children_assignment_combinations = 0
        # for correct_partial_result in inference_result:
        #     correct_children_assignment_combinations += 1
        #     (_, valid_candidates) = correct_partial_result
        #     total_candidates_found += len(valid_candidates)
        self.log(f"Total candidates found: {total_candidates_found}")
        # self.log(f"Correct children combinations found: {correct_children_assignment_combinations}")

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
        # TODO: Initialize the overlapping vertices correctly even if the space is passed manually
        self.overlapping_vertices = set()
        if initial_mapping is None:
            if current_initial_matching_mode == InitialMatchingMode.PLOID:
                initial_mapping = {x: [x] for x in processed_graph.probands}
            elif current_initial_matching_mode == InitialMatchingMode.INDIVIDUAL:
                initial_mapping = {x: [2 * (x // 2), (2 * (x // 2) + 1)] for x in processed_graph.probands}
                individual_mentioned = defaultdict(int)
                for vertex in initial_mapping:
                    vertex_individual = vertex // 2
                    individual_mentioned[vertex_individual] += 1
                self.overlapping_vertices = {x for x in individual_mentioned if individual_mentioned[x] == 2}
            else:
                raise Exception("Unsupported initial matching mode")
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
            coalescent_tree_vertex_to_subtrees[vertex] = self.get_subtree_matcher_for_coalescent_tree_proband(
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
            root_vertices = self.coalescent_tree.get_roots_for_clade(clade)
            for root_vertex in root_vertices:
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
                        print(f"Time spent on validation: {time_taken}. "
                              f"Average time taken: {time_taken / len(alignments)}")
                result_dict[root_vertex] = clade_alignments
        end_time = time.time()
        if print_enabled:
            print(f"Time spent on building and filtering the results: {end_time - start_time}")
        return result_dict

    def verify_valid_alignment(self, alignment: dict, root_vertex=None):
        """!
        @brief Verifies that the given alignment is correct
        """
        if len(alignment) != len(set(alignment.values())):
            return False
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
            self.add_edges_to_mrca_from_descendants(network_graph, parent_pedigree, children_pedigree)
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
                                        children_assignments=[{child: y}]) for x, y in vertex_subtree_dict[child].items()}
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
        if print_enabled:
            print(f"Inference time: {time_taken}")
            print("####################")
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
            candidate_subtree_matcher_dictionary[focal_vertex_candidate] = SubtreeMatcher(root_coalescent_tree=focal_vertex,
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
        if print_enabled:
            print(f"Building the result {result_build_end - result_build_start}")
            print(f"There are {len(candidate_subtree_matcher_dictionary)} resulting assignments")
        return candidate_subtree_matcher_dictionary

    def get_subtree_matcher_for_coalescent_tree_proband(self, proband: int, proband_pedigree_ids: [int]):
        """
        @brief Helper function returning the valid sub-alignments for a proband (which is simply one identity alignment).
        @param proband The id of the proband in the coalescent tree
        @param proband_pedigree_ids The ids of the corresponding vertices in the pedigree
        """
        return {proband_pedigree_id: SubtreeMatcher(root_coalescent_tree=proband, root_pedigree=proband_pedigree_id)
                for proband_pedigree_id in proband_pedigree_ids}

    # ----------------------------------------- Alignment logic ------------------------------------------------------

    def verify_pmrca_for_vertex_pair(self, first_vertex: int, second_vertex: int, common_ancestor: int):
        """!
        @brief This function verifies that the given vertex is a potential mrca for the given pair of vertices.
        @param first_vertex The first vertex from the pair.
        @param second_vertex The second vertex from the pair.
        @param common_ancestor A common ancestor of the given pair of vertices that is being verified.
        @return True if the given vertex is potential mrca for the given pair of vertices, False otherwise.
        """
        return not (len(self.pedigree.vertex_to_ancestor_map[first_vertex][common_ancestor]) == 1 and
                    self.pedigree.vertex_to_ancestor_map[first_vertex][common_ancestor] ==
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
        second_vertex_access = self.pedigree.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
        current_access = set(first_vertex_access).union(second_vertex_access)
        if len(current_access) > 2:
            return True
        third_vertex_access = self.pedigree.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
        for access_vertex in third_vertex_access:
            if access_vertex not in current_access:
                return True
        return False

    def extend_pair_to_triple(self, first_vertex: int, second_vertex: int, third_vertex: int,
                              potential_mrca_candidate: int):
        """!
        @brief This function verifies that the Hall condition holds for all the sets involving the third vertex.
        Before verifying the actual constraint, this function tries to check for the most common scenarios and avoid
        unnecessary set calculations.
        @param first_vertex The first vertex from the triplet.
        @param second_vertex The second vertex from the triplet.
        @param third_vertex The third vertex from the triplet.
        @param potential_mrca_candidate A common ancestor of the given triplet of vertices that satisfies the
        Hall's condition for the (first_vertex, second_vertex) pair.
        @return True if potential mrca candidate is indeed a potential mrca, False otherwise.
        """
        first_vertex_access = self.pedigree.vertex_to_ancestor_map[first_vertex][potential_mrca_candidate]
        second_vertex_access = self.pedigree.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
        third_vertex_access = self.pedigree.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
        if first_vertex_access == second_vertex_access == third_vertex_access:
            return False
        if len(third_vertex_access) == 1:
            if third_vertex_access == first_vertex_access:
                return False
            if third_vertex_access == second_vertex_access:
                return False
        return self.triplet_condition_holds(first_vertex, second_vertex, third_vertex, potential_mrca_candidate)

    def verify_additional_constraints(self, partially_assigned_vertices: [int], next_vertex: int, potential_mrca: int):
        """!
        @brief
        @param partially_assigned_vertices:
        @param next_vertex:
        @param potential_mrca:
        @return
        """
        if len(partially_assigned_vertices) == 2:
            [first_vertex, second_vertex] = partially_assigned_vertices
            return self.extend_pair_to_triple(first_vertex, second_vertex, next_vertex, potential_mrca)
        raise Exception("Unimplemented")

    def apply_additional_constraints_for_partial_result(self, coalescent_vertex_to_candidates: {int: [int]},
                                                        partially_assigned_vertices: [int], potential_mrcas: [int],
                                                        next_vertex: int,
                                                        ):
        """!
        @brief
        @param coalescent_vertex_to_candidates:
        @param partially_assigned_vertices:
        @param potential_mrcas:
        @param next_vertex:
        @return:
        """
        next_vertex_candidates = coalescent_vertex_to_candidates[next_vertex]
        next_vertex_result = []
        for candidate in next_vertex_candidates:
            candidate_ancestors = set(self.pedigree.vertex_to_ancestor_map[candidate])
            candidate_result = []
            for potential_mrca in potential_mrcas:
                if (potential_mrca not in candidate_ancestors or
                        not self.verify_additional_constraints(partially_assigned_vertices, next_vertex,
                                                               potential_mrca)):
                    continue
                candidate_result.append(potential_mrca)
            if candidate_result:
                next_partially_assigned_vertices = list(partially_assigned_vertices)
                next_partially_assigned_vertices.append(candidate)
                next_vertex_result.append((next_partially_assigned_vertices, candidate_result))
        return next_vertex_result

    def get_single_vertex_subsets(self, vertex: int, vertex_ancestors):
        """!
        @brief A helper function returning the partial sub-assignment result for the specified vertex with the
        given ancestors.
        @param vertex
        @param vertex_ancestors
        @return The sub-assignment tuple for a single vertex used by the alignment algorithm
        """
        return [(ancestor, [(1, frozenset(self.pedigree.vertex_to_ancestor_map[vertex][ancestor]))]) for
                ancestor in vertex_ancestors]

    def get_single_vertex_verified_subset_tuple(self, vertex: int, vertex_ancestors):
        """!
        @brief
        @param vertex:
        @param vertex_ancestors:
        @return:
        """
        return (vertex,), self.get_single_vertex_subsets(vertex, vertex_ancestors)

    def filter_common_ancestors_for_vertex_pair(self, first_candidate, second_candidate):
        """!
        @brief
        @param first_candidate:
        @param second_candidate:
        @return:
        """
        verified_ancestors = []
        if first_candidate not in self.pedigree.parents_map or second_candidate not in self.pedigree.parents_map:
            return verified_ancestors
        while self.pedigree.parents_map[first_candidate] == self.pedigree.parents_map[second_candidate]:
            verified_ancestors.extend(self.pedigree.parents_map[first_candidate])
            [first_candidate, second_candidate] = self.pedigree.parents_map[first_candidate]
            if first_candidate not in self.pedigree.parents_map or second_candidate not in self.pedigree.parents_map:
                return verified_ancestors
        first_ancestors = self.pedigree.vertex_to_ancestor_map[first_candidate]
        second_ancestors = self.pedigree.vertex_to_ancestor_map[second_candidate]
        if len(second_ancestors) > len(first_ancestors):
            first_ancestors, second_ancestors = second_ancestors, first_ancestors
        for ancestor in first_ancestors:
            if (ancestor in second_ancestors and
                    self.verify_pmrca_for_vertex_pair(first_candidate, second_candidate,
                                                      ancestor)):
                verified_ancestors.append(ancestor)
        return verified_ancestors

    def get_pmracs_for_vertex_pair(self, first: int, second: int, coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first:
        @param second:
        @param coalescent_vertex_to_candidates:
        @return
        """
        result = defaultdict(list)
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        if len(first_vertex_candidates) > len(second_vertex_candidates):
            first_vertex_candidates, second_vertex_candidates = second_vertex_candidates, first_vertex_candidates
        for first_vertex_candidate in first_vertex_candidates:
            for second_vertex_candidate in second_vertex_candidates:
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(
                    first_vertex_candidate,
                    second_vertex_candidate)
                for verified_ancestor in verified_ancestors:
                    result[verified_ancestor].append((first_vertex_candidate, second_vertex_candidate))
        return result

    def get_pmracs_for_vertex_triple_dfs(self, first: int, second: int, third: int,
                                         coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first:
        @param second:
        @param third:
        @param coalescent_vertex_to_candidates:
        @return:
        """
        result = []
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        third_vertex_candidates = coalescent_vertex_to_candidates[third]
        for first_vertex_candidate in first_vertex_candidates:
            first_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[first_vertex_candidate]
            if not first_vertex_candidate:
                continue
            for second_vertex_candidate in second_vertex_candidates:
                second_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[second_vertex_candidate]
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(first_vertex_candidate,
                                                                                  second_vertex_candidate)
                if not verified_ancestors:
                    continue
                for third_vertex_candidate in third_vertex_candidates:
                    third_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[third_vertex_candidate]
                    if not third_vertex_ancestors:
                        continue
                    if (self.pedigree.parents_map[first_vertex_candidate] ==
                            self.pedigree.parents_map[second_vertex_candidate] ==
                            self.pedigree.parents_map[third_vertex_candidate]):
                        fully_verified_ancestors = self.pedigree.parents_map[first_vertex_candidate]
                    else:
                        fully_verified_ancestors = []
                        for verified_ancestor in verified_ancestors:
                            if (verified_ancestor in third_vertex_ancestors and
                                    self.extend_pair_to_triple(first_vertex_candidate, second_vertex_candidate,
                                                               third_vertex_candidate, verified_ancestor)):
                                fully_verified_ancestors.append(verified_ancestor)
                    if fully_verified_ancestors:
                        result.append(((first_vertex_candidate, second_vertex_candidate, third_vertex_candidate),
                                       fully_verified_ancestors
                                       ))
        return result

    def get_pmracs_for_vertex_triple_iterative(self, first: int, second: int, third: int,
                                               coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first: The first coalescent vertex
        @param second:The second coalescent vertex
        @param third: The third coalescent vertex
        @param coalescent_vertex_to_candidates: The map returning the candidate list for a given coalescent vertex
        @return:
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
        shared_common_ancestors = set(first_third_pair_result.keys()).intersection(
            first_third_pair_result.keys())
        result = defaultdict(list)
        for shared_common_ancestor in shared_common_ancestors:
            shared_first_vertex_assignments = (frozenset(first_second_dict[shared_common_ancestor]).
                                               intersection(first_third_dict[shared_common_ancestor]))
            for first_candidate in shared_first_vertex_assignments:
                for second_candidate in first_second_dict[shared_common_ancestor][first_candidate]:
                    for third_candidate in first_third_dict[shared_common_ancestor][first_candidate]:
                        triplet_tuple = get_triplet_tuple(first_candidate, second_candidate, third_candidate)
                        # TODO: Add common parents speed-up
                        if (self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, shared_common_ancestor)
                                and self.triplet_condition_holds(first_candidate, second_candidate,
                                                                 third_candidate, shared_common_ancestor)):
                            result[shared_common_ancestor].append(triplet_tuple)

        # result = []
        # for (first_candidate, first_second_partial_results) in first_second_dict.items():
        #     if first_candidate not in first_third_dict:
        #         continue
        #     for ((_, second_candidate), verified_ancestors_second) in first_second_partial_results:
        #         assert first_candidate == _
        #         verified_ancestors_second_parents_restricted = None
        #         second_candidate_ancestors = self.pedigree.vertex_to_ancestor_map[second_candidate]
        #         for first_third_partial_result in first_third_dict[first_candidate]:
        #             ((__, third_candidate), verified_ancestors_third) = first_third_partial_result
        #             assert first_candidate == __
        #             resulting_ancestors = verified_ancestors_second
        #             # If the second and the third candidates happen to have the same parents, we can restrict
        #             # the candidates space that should be searched. Specifically, the resulting pmracs
        #             # must belong to the intersection of the parents' ancestors sets or be the parents themselves
        #             if self.pedigree.parents_map[second_candidate] == self.pedigree.parents_map[third_candidate]:
        #                 if self.pedigree.parents_map[first_candidate] == self.pedigree.parents_map[second_candidate]:
        #                     result.append(((first_candidate, second_candidate, third_candidate),
        #                                    self.pedigree.parents_map[first_candidate]))
        #                     continue
        #                 if verified_ancestors_second_parents_restricted is None:
        #                     [first_parent, second_parent] = self.pedigree.parents_map[second_candidate]
        #                     parents_set = {first_parent, second_parent}
        #                     first_parent_ancestors = self.pedigree.vertex_to_ancestor_map[first_parent]
        #                     second_parent_ancestors = self.pedigree.vertex_to_ancestor_map[second_parent]
        #                     parents_common_ancestors = set(first_parent_ancestors).intersection(
        #                         set(second_parent_ancestors))
        #                     verified_ancestors_second_parents_restricted = (
        #                         (parents_set.union(parents_common_ancestors)).intersection(verified_ancestors_second))
        #                 resulting_ancestors = verified_ancestors_second_parents_restricted
        #             # The resulting pmracs should belong to the both partial pmracs sets
        #             resulting_ancestors = [x for x in verified_ancestors_third
        #                                    if x in resulting_ancestors and
        #                                    self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, x) and
        #                                    self.triplet_condition_holds(first_candidate, second_candidate,
        #                                                                 third_candidate, x)]
        #             if resulting_ancestors:
        #                 result.append(((first_candidate, second_candidate, third_candidate), resulting_ancestors))
        return result

    def get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
                                coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief This function calculates the potential most recent common ancestor (pmrca) for the given vertices.
        TODO: Give the definition of a pmrca taking into account various optimizations that have been implemented
        @param vertices_coalescent_ids The ids for the coalescent vertices.
        @param coalescent_vertex_to_candidates: A dictionary mapping an id of a coalescent vertex to the list of
        ids of pedigree vertices that can be assigned to this coalescent vertex.
        @return All the valid pmracs for the given vertices.
        """
        vertices_length = len(vertices_coalescent_ids)
        vertices_coalescent_ids = sorted(vertices_coalescent_ids,
                                         key=lambda child: len(
                                             coalescent_vertex_to_candidates[child]), reverse=False)
        if vertices_length == 2:
            result = self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                     coalescent_vertex_to_candidates)
        elif vertices_length == 3:
            result = self.get_pmracs_for_vertex_triple_iterative(vertices_coalescent_ids[0],
                                                                 vertices_coalescent_ids[1],
                                                                 vertices_coalescent_ids[2],
                                                                 coalescent_vertex_to_candidates)
        # elif vertices_length < 5:
        #     result = self.get_pmracs_for_vertices_dfs(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        else:
            separate_results = []
            start_index = 0
            if len(vertices_coalescent_ids) % 2 == 1:
                separate_results.append(self.get_pmracs_for_vertex_triple_iterative(
                    vertices_coalescent_ids[0],
                    vertices_coalescent_ids[1],
                    vertices_coalescent_ids[2],
                    coalescent_vertex_to_candidates
                ))
                start_index = 3
            separate_results.extend([self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[i],
                                                                     vertices_coalescent_ids[i + 1],
                                                                     coalescent_vertex_to_candidates)
                                     for i in range(start_index, len(vertices_coalescent_ids), 2)])
            common_keys = frozenset.intersection(*(frozenset(d.keys()) for d in separate_results))
            separate_results = [ {key: d[key] for key in common_keys} for d in separate_results]
            result = separate_results[0]
            for next_result in separate_results[1:]:
                new_result = dict()
                for ancestor_candidate, children_assignments in result.items():
                    extended_children_assignments = []
                    for children_assignment in children_assignments:
                        for next_children_assignment in next_result[ancestor_candidate]:
                            joined_children_assignment = children_assignment + next_children_assignment
                            if (len(frozenset(next_children_assignment).union(children_assignment))
                                    != len(children_assignment) + len(next_children_assignment)):
                                continue
                            if (len(joined_children_assignment) > 3 or
                                    self.verify_pmrca_for_vertices(ancestor_candidate, joined_children_assignment)):
                                extended_children_assignments.append(joined_children_assignment)
                    if extended_children_assignments:
                        new_result[ancestor_candidate] = extended_children_assignments
                result = new_result
        # result = self.get_pmrcas_without_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        # result = self.get_pmracs_for_vertices_with_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        # result = self.get_mrcas_without_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        # for assignment in result:
        #     (assigned_children, verified_candidates_list) = assignment
        #     self.inference_cache[assigned_children] = verified_candidates_list
        if not result:
            raise Exception("Failed")
        return vertices_coalescent_ids, result

    def get_pmrcas_without_memory(self, vertices_coalescent_ids: [int], coalescent_vertex_to_candidates: {int: [int]}):
        vertices_length = len(vertices_coalescent_ids)
        first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
        partial_result = []
        for first_candidate in first_candidates:
            # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
            # in the coalescent tree for which we are making the inference
            first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(first_candidate)
                                         if len(self.pedigree.children_map[x]) >= vertices_length]
            if first_candidate_ancestors:
                partial_result.append(([first_candidate], first_candidate_ancestors))
        for next_coalescent_id in vertices_coalescent_ids[1:]:
            next_result = []
            for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
                next_candidate_ancestors = self.pedigree.get_vertex_ancestors(next_candidate)
                for verified_children_partial_assignment in partial_result:
                    (assigned_children, verified_candidates_list) = verified_children_partial_assignment
                    new_assigned_children_list = assigned_children + [next_candidate]
                    new_verified_candidates_list = []
                    for candidate in verified_candidates_list:
                        if candidate not in next_candidate_ancestors:
                            continue
                        if self.verify_pmrca_for_vertices(candidate, new_assigned_children_list):
                            new_verified_candidates_list.append(candidate)
                    if new_verified_candidates_list:
                        next_result.append((new_assigned_children_list, new_verified_candidates_list))
            partial_result = next_result
        return partial_result

    def get_pmracs_for_vertices_with_maximum_flow_iterative(self, vertices_coalescent_ids: [int],
                                                            coalescent_vertex_to_candidates: {int: [int]}):
        vertices_length = len(vertices_coalescent_ids)
        first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
        partial_result = []
        for first_candidate in first_candidates:
            # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
            # in the coalescent tree for which we are making the inference
            first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(first_candidate)
                                         if len(self.pedigree.children_map[x]) >= vertices_length]
            if first_candidate_ancestors:
                partial_result.append((first_candidate, first_candidate_ancestors))
        for next_coalescent_id in vertices_coalescent_ids[1:]:
            next_result = []
            if print_enabled:
                print(f"Partial result: {len(partial_result)}")
                print(f"Candidates for the next vertex: {len(coalescent_vertex_to_candidates[next_coalescent_id])}")
                print(
                    f"The whole search space is {len(partial_result) * len(coalescent_vertex_to_candidates[next_coalescent_id])}")
            for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
                next_candidate_ancestors = self.pedigree.get_vertex_ancestors(next_candidate)
                for verified_children_partial_assignment in partial_result:
                    (assigned_children, verified_candidates_list) = verified_children_partial_assignment
                    next_verified_candidates_list = []
                    extended_assigned_children = assigned_children + (next_candidate,)
                    for verified_candidate in verified_candidates_list:
                        if verified_candidate not in next_candidate_ancestors:
                            continue
                        if self.verify_mrca_for_vertices(verified_candidate, extended_assigned_children):
                            next_verified_candidates_list.append(verified_candidate)
                    if next_verified_candidates_list:
                        next_result.append((extended_assigned_children, next_verified_candidates_list))
            if print_enabled:
                print(f"The resulting number of correct assignments is: {len(next_result)}")
            partial_result = next_result
        if not partial_result:
            raise Exception("No valid assignments found")
        return partial_result

    def get_pmracs_for_vertices_dfs(self, vertices_coalescent_ids: [int],
                                    coalescent_vertex_to_candidates: {int: [int]}):
        result = []
        vertices_length = len(vertices_coalescent_ids)
        current_index = 0
        current_partial_assignment = dict()
        partial_result = dict()
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
        return result

    # Returns all the potential common ancestors for the specified vertices
    # Note that not all the potential common ancestors are MRCAs
    def get_pmracs_for_vertices_with_memory(self, vertices_coalescent_ids: [int],
                                            coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief Calculates the pmracs for the given vertices in an iterative manner storing the calculations done
        during the verification of the previous sub-assignment.
        """
        vertices_length = len(vertices_coalescent_ids)
        # partial_result elements are lists whose elements
        # have the format: (partial_assignment: tuple, [(candidate, [(preimage_length, image)])])
        first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
        partial_result = Cache(size_limit=2 * 1024 * 1000000, eviction_policy="FIFO")
        for first_candidate in first_candidates:
            # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
            # in the coalescent tree for which we are making the inference
            first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(first_candidate)
                                         if len(self.pedigree.children_map[x]) >= vertices_length]
            if first_candidate_ancestors:
                (assigned_children, verified_candidates_list) = self.get_single_vertex_verified_subset_tuple(
                    first_candidate,
                    first_candidate_ancestors)
                partial_result[assigned_children] = verified_candidates_list
                # partial_result.append(self.get_single_vertex_verified_subset_tuple(first_candidate,
                #                                                                  first_candidate_ancestors))
        for next_coalescent_id in vertices_coalescent_ids[1:]:
            next_result = Cache(size_limit=2 * 1024 * 1000000, eviction_policy="FIFO")
            if print_enabled:
                print(f"Partial result: {len(partial_result)}")
                print(f"Candidates for the next vertex: {len(coalescent_vertex_to_candidates[next_coalescent_id])}")
                print(
                    f"The whole search space is {len(partial_result) * len(coalescent_vertex_to_candidates[next_coalescent_id])}")
            for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
                next_candidate_ancestors = self.pedigree.get_vertex_ancestors(next_candidate)
                for assigned_children in partial_result.iterkeys():
                    verified_candidates_list = partial_result[assigned_children]
                    next_verified_candidates_list = []
                    extended_assigned_children = assigned_children + (next_candidate,)
                    for verified_candidate_tuple in verified_candidates_list:
                        (verified_candidate, verified_subsets) = verified_candidate_tuple
                        if verified_candidate not in next_candidate_ancestors:
                            continue
                        is_correct_assignment = True
                        verified_candidate_access = self.pedigree.vertex_to_ancestor_map[next_candidate][
                            verified_candidate]
                        new_subsets = []
                        for verified_subset_tuple in verified_subsets:
                            (preimage_size, image) = verified_subset_tuple
                            image: frozenset
                            new_image = image.union(verified_candidate_access)
                            if len(new_image) < preimage_size + 1:
                                is_correct_assignment = False
                                break
                            new_subsets.append((preimage_size + 1, new_image))
                        if not is_correct_assignment:
                            continue
                        new_subsets.extend(verified_subsets)
                        new_subsets.append(
                            (1, set(self.pedigree.vertex_to_ancestor_map[next_candidate][verified_candidate])))
                        next_verified_candidates_list.append((verified_candidate, new_subsets))
                        # next_result.append((extended_assigned_children, new_subsets))
                    if next_verified_candidates_list:
                        next_result[extended_assigned_children] = next_verified_candidates_list
            if print_enabled:
                print(f"The resulting number of correct assignments is: {len(next_result)}")
            partial_result = next_result
        # Refactored code
        # alignment_result = [(key, partial_result[key]) for key in partial_result]
        result = [(assigned_children, [candidate[0] for candidate in partial_result[assigned_children]]) for
                  assigned_children in partial_result]
        # result = []
        # for assignment in alignment_result:
        #     (assigned_children, verified_candidates_list) = assignment
        #     candidates = [verified_candidate_tuple[0] for verified_candidate_tuple in verified_candidates_list]
        #     result.append((assigned_children, candidates))
        if not result:
            raise Exception("No valid assignments found")
        return result

    def verify_mrca_for_vertices(self, mrca: int, descendant_vertices: [int]):
        # A speed-up for the case with two vertices
        if len(descendant_vertices) == 2 and \
                self.pedigree.vertex_to_ancestor_map[descendant_vertices[0]][mrca] != \
                self.pedigree.vertex_to_ancestor_map[descendant_vertices[1]][mrca]:
            return True
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        for descendant in descendant_vertices:
            network_graph.add_edge(source_vertex_label, descendant, capacity=1)
            self.add_edges_to_mrca_from_descendants(network_graph, mrca, descendant_vertices)
        return networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=mrca) == len(
            descendant_vertices)

    def verify_pmrca_for_vertices(self, pmrca: int, descendant_vertices: [int]):
        # A speed-up for the case with two vertices
        if len(descendant_vertices) == 2 and \
                self.pedigree.vertex_to_ancestor_map[descendant_vertices[0]][pmrca] != \
                self.pedigree.vertex_to_ancestor_map[descendant_vertices[1]][pmrca]:
            return True
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        for descendant in descendant_vertices:
            network_graph.add_edge(source_vertex_label, descendant, capacity=1)
            for access_vertex in self.pedigree.vertex_to_ancestor_map[descendant][pmrca]:
                network_graph.add_edge(descendant, access_vertex, capacity=1)
                network_graph.add_edge(access_vertex, pmrca, capacity=1)
        return networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=pmrca) == len(
            descendant_vertices)

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
