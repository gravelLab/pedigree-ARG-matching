"""!
@file subtree_matcher.py
@brief Brief description here
"""

from __future__ import annotations
from itertools import product


class SubtreeMatcher:
    """
    This class represents multiple matchings between a subtree of a coalescent tree in the ARG and the pedigree.
    # TODO: Add some documentation about the expected behaviour

    Attributes:
        root_coalescent_tree (int): Represents the vertex id of the root of the coalescent tree
        root_pedigree (int): Represents the vertex in the pedigree to which root_coalescent_tree
        can be mapped.
        children_assignments ([{int: SubtreeMatcher}]): The dictionaries representing the assignments to the children
        vertices which make root_pedigree a valid assignment to root_coalescent_tree
    """

    def __init__(self, root_coalescent_tree: int, root_pedigree: int,
                 children_assignments: [{int: SubtreeMatcher}] = None,
                 overlapping_ploids_assignments: [{int: [int]}] = None):
        self.root_coalescent_tree = root_coalescent_tree
        self.root_pedigree = root_pedigree
        self.children_assignments = children_assignments
        self.overlapping_ploids_assignments = overlapping_ploids_assignments
        self.subtree_alignments = None

    def get_all_subtree_alignments(self):
        results = []
        root_assignment_dict = {self.root_coalescent_tree: self.root_pedigree}
        # If the given vertex has no children, the resulting alignment is a dictionary with one key-value pair
        if self.children_assignments is None:
            return [root_assignment_dict]
        # If there are children assignments, loop over all the corresponding children alignments
        for children_assignment in self.children_assignments:
            children_alignments = [x.get_all_subtree_alignments() for x in children_assignment.values()]
            for children_dictionaries in product(*children_alignments):
                new_result = dict(root_assignment_dict)
                for dictionary in children_dictionaries:
                    new_result.update(dictionary)
                results.append(new_result)
        self.subtree_alignments = results
        return results


class AlignmentsIterator:
    def __init__(self, root_coalescent_tree: int, vertex_subtree_matchers: {int: [SubtreeMatcher]},
                 vertex_to_dictionaries: {int: [dict]} = None):
        if root_coalescent_tree not in vertex_subtree_matchers:
            raise Exception("The root coalescent tree has no available matchings")
        self.root_coalescent_tree = root_coalescent_tree
        self.root_matchers = vertex_subtree_matchers[root_coalescent_tree]
        self.vertex_subtree_matchers = vertex_subtree_matchers
        self.subtree_matcher_counter = 0
        self.subtree_matchers_iterators = None
        # The cache that stores all the alignments for the given vertex in the coalescent tree
        if vertex_to_dictionaries is None:
            vertex_to_dictionaries = dict()
        self.vertex_to_dictionaries = vertex_to_dictionaries
