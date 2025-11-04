from __future__ import annotations
from itertools import product


class SubtreeMatcher:
    """
    This class represents multiple matchings between a subtree of a coalescent tree in the ARG and the pedigree.

    Attributes:
        root_coalescent_tree (int): Represents the vertex id of the root of the coalescent tree
        root_pedigree (int): Represents the vertex in the pedigree to which root_coalescent_tree
        can be mapped.
        children_assignments ([{int: SubtreeMatcher}]): The dictionaries representing the assignments to the children
        vertices which make root_pedigree a valid assignment to root_coalescent_tree
    """

    def __init__(self, root_coalescent_tree: int, root_pedigree: int,
                 children_assignments: [{int: SubtreeMatcher}] = None):
        """
        Initializes the SubtreeMatcher object.

        Args:

            root_coalescent_tree (int): Represents the vertex id of the root of the coalescent tree
            root_pedigree (int): Represents the vertex in the pedigree to which root_coalescent_tree
            can be mapped.
            children_assignments ([{int: SubtreeMatcher}]): The dictionaries representing the assignments
            to the children vertices which make root_pedigree a valid assignment to root_coalescent_tree
        """
        self.root_coalescent_tree = root_coalescent_tree
        self.root_pedigree = root_pedigree
        self.children_assignments = children_assignments
        self.subtree_alignments = None

    def get_all_subtree_alignments(self):
        """
        Builds all the distinct valid alignment for the coalescent vertex subclade.
        """
        root_assignment_dict = {self.root_coalescent_tree: self.root_pedigree}

        if self.children_assignments is None:
            yield root_assignment_dict
            return

        for children_assignment in self.children_assignments:
            children_alignments = [x.get_all_subtree_alignments() for x in children_assignment.values()]
            for children_dicts in product(*children_alignments):
                new_result = dict(root_assignment_dict)
                for d in children_dicts:
                    new_result.update(d)
                yield new_result
