"""!
@file subtree_matcher.py
@brief Brief description here
"""


class SubtreeMatcher:
    """
    This class represents a matching between a subtree of a coalescent tree in the ARG and the pedigree.

    Attributes:
        root_coalescent_tree (int): Represents the vertex id of the root of the coalescent tree
        root_pedigree (int): Represents the vertex if of the pedigree vertex to which root_coalescent_tree
         is mapped.
        subtree_matchers (dict): A dictionary representing the assignments to the children vertices which make
        root_pedigree a valid assignment to root_coalescent_tree
    """

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
