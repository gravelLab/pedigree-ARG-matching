
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from genealogical_graph import CoalescentTree

import time
import datetime

logs_enabled = False


class MatcherLogger:

    def __init__(self):
        if logs_enabled:
            time_in_milliseconds = datetime.datetime.now().timestamp() * 1000
            log_filename = f"logs/log_{int(time_in_milliseconds)}.txt"
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

    def close(self):
        if logs_enabled:
            self.file.close()


logger = MatcherLogger()


class SubtreeMatcher:

    def __init__(self, root_coalescent_tree: int, root_pedigree: int,
                 subtrees_matchers: dict = None):
        self.root_coalescent_tree = root_coalescent_tree
        self.root_pedigree = root_pedigree
        self.subtree_matchers = subtrees_matchers
        self.dict_assignment = None

    def get_dict_assignment(self):
        if self.subtree_matchers is None:
            return {self.root_pedigree: self.root_coalescent_tree}
        if self.dict_assignment is None:
            # TODO: Think whether we need to consider correctness of the assignment here
            children_assignments = [x.get_dict_assignment() for x in self.subtree_matchers]
            merged_dictionary = dict()
            for dictionary in children_assignments:
                merged_dictionary.update(dictionary)
            self.dict_assignment = merged_dictionary
        return self.dict_assignment


def get_subtrees_from_children(focal_vertex: int, coalescent_tree: CoalescentTree,
                               pedigree: PotentialMrcaProcessedGraph,
                               vertex_subtree_dict: {int: [SubtreeMatcher]}):
    focal_vertex_children: [int] = coalescent_tree.children_map[focal_vertex]
    if len(focal_vertex_children) == 0:
        raise Exception("Isolated vertex in the coalescent tree")
    if len(focal_vertex_children) == 1:
        # This situation must never happen with real coalescent trees. We could have just 'squashed' the paths
        # consisting of vertices having just one child, but we will lose some information this way.
        # This usually leads to a much bigger search space that we wouldn't have gotten with normal coalescent tree
        # Therefore, we will cheat here
        [child] = focal_vertex_children
        child_potential_pedigree_values = {x.root_pedigree for x in vertex_subtree_dict[child]}
        result = [SubtreeMatcher(root_coalescent_tree=focal_vertex, root_pedigree=focal_vertex,
                                 subtrees_matchers={child: x})
                  for x in child_potential_pedigree_values]
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
    inference_result = pedigree.get_pmracs_for_vertices(
        coalescent_vertex_to_candidates=child_candidate_subtree_matchers_matrix,
        vertices_coalescent_ids=focal_vertex_children
    )
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
