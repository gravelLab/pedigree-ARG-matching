from path_matching.constraint import Constraint
from path_matching.variable import Variable


class CommonAncestorConstraint(Constraint):

    def __init__(self, common_ancestor_coalescent_tree_id: int, variables: [Variable]):
        self.common_ancestor = common_ancestor_coalescent_tree_id
        self.variables = variables

    def verify(self):
        result = False
        ancestors_candidate_list = [variable.ancestor_pedigree_candidates_map.keys() for variable in self.variables]
        common_candidates = set.intersection(*ancestors_candidate_list)
        for variable in self.variables:
            narrowed_candidates_map = {key: value for key, value in variable.ancestor_pedigree_candidates_map.items() if key in common_candidates}
            result = result or len(narrowed_candidates_map) < len(variable.ancestor_pedigree_candidates_map)
            variable.ancestor_pedigree_candidates_map = narrowed_candidates_map
        return result
