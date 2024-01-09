import datetime
import os
import time

from divide_and_conquer.subtree_matcher import *

logs_enabled = True
logs_default_directory_name = "logs"

print_enabled = True


class MatcherLogger:

    def __init__(self, logs_directory_name=logs_default_directory_name):
        if logs_enabled:
            time_in_milliseconds = datetime.datetime.now().timestamp() * 1000
            log_filename = f"{logs_directory_name}/log_{int(time_in_milliseconds)}.txt"
            if not os.path.exists(logs_directory_name):
                os.makedirs(logs_directory_name)
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

    def __init__(self, processed_graph: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree,
                 logger: MatcherLogger, initial_mapping: dict = None):
        self.pedigree = processed_graph
        self.coalescent_tree = coalescent_tree
        if initial_mapping is None:
            initial_mapping = {x: x for x in processed_graph.probands}
        self.initial_mapping = initial_mapping
        self.logger = logger

    def find_mapping(self):
        coalescent_tree_vertex_to_subtrees = dict()
        for vertex in self.coalescent_tree.levels[0]:
            coalescent_tree_vertex_to_subtrees[vertex] = get_subtree_matcher_for_coalescent_tree_proband(
                vertex, self.initial_mapping[vertex])
        for level in self.coalescent_tree.levels[1:]:
            for vertex in level:
                vertex_subtrees = self.get_subtrees_from_children(vertex,
                                                                  coalescent_tree_vertex_to_subtrees)
                coalescent_tree_vertex_to_subtrees[vertex] = vertex_subtrees
        # root_assignments = [coalescent_tree_vertex_to_subtrees[x] for x in self.coalescent_tree.levels[-1]]
        return coalescent_tree_vertex_to_subtrees

    def get_subtrees_from_children(self, focal_vertex: int,
                                   vertex_subtree_dict: {int: [SubtreeMatcher]}):
        focal_vertex_children: [int] = self.coalescent_tree.children_map[focal_vertex]
        if len(focal_vertex_children) == 0:
            raise Exception("Isolated vertex in the coalescent tree")
        if len(focal_vertex_children) == 1:
            # This situation must never happen with real coalescent trees. We could have just 'squashed' the paths
            # consisting of vertices having just one child, but we will lose some information this way.
            # This usually leads to a much bigger search space that we wouldn't have gotten with normal coalescent tree
            # Therefore, we will cheat here
            [child] = focal_vertex_children
            # The correct one
            result = [SubtreeMatcher(root_coalescent_tree=focal_vertex, root_pedigree=x.root_pedigree,
                                     subtrees_matchers={child: x.root_pedigree}) for x in vertex_subtree_dict[child]]
            return result
        child_candidate_subtree_matchers_matrix = {child: dict() for child in focal_vertex_children}
        child_candidate_subtree_matchers_matrix: {int: {int: [SubtreeMatcher]}}
        for child in focal_vertex_children:
            for subtree_matcher in vertex_subtree_dict[child]:
                subtree_matcher: SubtreeMatcher
                child_candidate = subtree_matcher.root_pedigree
                if child_candidate not in child_candidate_subtree_matchers_matrix[child]:
                    child_candidate_subtree_matchers_matrix[child][child_candidate] = list()
                child_candidate_subtree_matchers_matrix[child][child_candidate].append(subtree_matcher)
        # print(f"Doing the inference for the vertex {focal_vertex} which has {len(focal_vertex_children)}")
        # for child in focal_vertex_children:
        #     print(f"{child}: {len(vertex_subtree_dict[child])}")
        if print_enabled:
            print("####################")
            print(f"Inference for {focal_vertex} ({self.pedigree.vertex_to_level_map[focal_vertex]}),"
                  f" there are {len(focal_vertex_children)} children")
            for child in focal_vertex_children:
                print(f"{child}: {len(child_candidate_subtree_matchers_matrix[child])}")
        inference_start = time.time()
        inference_result = self.pedigree.get_pmracs_for_vertices(
            coalescent_vertex_to_candidates=child_candidate_subtree_matchers_matrix,
            vertices_coalescent_ids=focal_vertex_children
        )
        if print_enabled:
            print("####################")
        inference_end = time.time()
        time_taken = inference_end - inference_start
        self.logger.log_vertex_inference_time(self.coalescent_tree, time_taken, focal_vertex, focal_vertex_children,
                                              child_candidate_subtree_matchers_matrix, inference_result)
        result = []
        for single_result in inference_result:
            (assigned_children, focal_vertex_candidates) = single_result
            assert len(assigned_children) == len(focal_vertex_children)
            children_dictionary = dict()
            for index in range(len(assigned_children)):
                children_dictionary[focal_vertex_children[index]] = assigned_children[index]
            result.extend([SubtreeMatcher(focal_vertex, x, children_dictionary) for x in focal_vertex_candidates])
        return result


def get_subtree_matcher_for_coalescent_tree_proband(proband: int, proband_pedigree_id: int):
    return [SubtreeMatcher(root_coalescent_tree=proband, root_pedigree=proband_pedigree_id)]
