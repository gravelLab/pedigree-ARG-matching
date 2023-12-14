from itertools import combinations

from path_matching.constraint import Constraint
from genealogical_graph import CoalescentTree
from path_matching.common_ancestor_constraint import CommonAncestorConstraint
from path_matching.path_processed_graph import PathProcessedGraph
from path_matching.variable import Variable


class PathGraphMatcher:
    class IntersectConstraint(Constraint):

        def __init__(self, coalescent_vertex_id: int, path_matcher):
            self.coalescent_vertex_id = coalescent_vertex_id
            self.coalescent_children = path_matcher.coalescent_tree.children_map[coalescent_vertex_id]
            self.path_matcher = path_matcher

        def verify(self):
            path_matcher: PathGraphMatcher
            result = False
            for child in self.coalescent_children:
                variables = self.path_matcher.find_common_ancestor_related_variables(child)
                child_candidate_edge_values_map = dict()
                for (a, b) in combinations(variables, 2):
                    candidate_map = self.path_matcher.path_processed_graph.path_common_ancestry_map[a][b]
                    for (lower_candidate, higher_candidate) in combinations(candidate_map.items(), 2):
                        if self.path_matcher.get_pedigree_vertex_level(
                                lower_candidate) >= self.path_matcher.get_pedigree_vertex_level(higher_candidate):
                            continue

            return result

    def set_trivial_initial_mapping(self):
        self.initial_mapping = {key: key for key in self.coalescent_tree.probands}

    def __init__(self, coalescent_tree: CoalescentTree,
                 path_processed_graph: PathProcessedGraph,
                 initial_mapping: dict = None):
        self.coalescent_tree = coalescent_tree
        self.path_processed_graph = path_processed_graph
        if initial_mapping is None:
            self.set_trivial_initial_mapping()
        else:
            self.initial_mapping = initial_mapping
        self.variables = dict()
        self.initialize_variables()
        self.constraints = set()
        self.initialize_constraints()

    def get_proband_ancestor_domain_map(self):
        result = dict()
        for proband_coalescent_tree in self.coalescent_tree.probands:
            proband_pedigree = self.initial_mapping[proband_coalescent_tree]
            proband_dict = dict()
            result[proband_pedigree] = proband_dict
            domain = [self.path_processed_graph.path_common_ancestry_map[proband_pedigree].items()]
        return result

    # def initialize_variables(self):
    #     ancestry_map = self.path_processed_graph.path_common_ancestry_map
    #     for level in self.coalescent_tree.levels:
    #         for vertex in level:
    #             coalescent_tree_children = self.coalescent_tree.children_map[vertex]
    #             child_to_domain_map = dict()
    #             for child in coalescent_tree_children:
    #                 if child in self.coalescent_tree.probands:
    #                     child_to_domain_map[child] = [child]
    #                 else:
    #                     pass
        # Processing the first level in the coalescent tree
        # for first_level_coalescent_tree_vertex in self.coalescent_tree.levels[1]:
        #     # Children in the pedigree
        #     children: [int] = [self.initial_mapping[x] for x in
        #                        self.coalescent_tree.children_map[first_level_coalescent_tree_vertex]]
        #     all_common_ancestors = self.path_processed_graph.get_common_ancestries_for_vertices(vertices=children)
        #     for child in children:
        #         child_dictionary = dict()
        #         for common_ancestor in all_common_ancestors:
        #             child_dictionary[common_ancestor] = self.path_processed_graph.path_map[common_ancestor][child]
        #             common_ancestor_child_variable = Variable(proband_pedigree=child,
        #                                                       ancestor_coalescent_tree=first_level_coalescent_tree_vertex,
        #                                                       domain=child_dictionary)
        #             self.variables[(first_level_coalescent_tree_vertex, child)] = common_ancestor_child_variable
        # # Processing the next level
        # for level in self.coalescent_tree.levels[2:]:
        #     for next_level_vertex in level:
        #         coalescent_tree_children = self.coalescent_tree.children_map[next_level_vertex]
        #         for child in coalescent_tree_children:

    def add_common_ancestor_constraint(self, unmapped_ancestor: int):
        unmapped_ancestor_proband_descendants = self.coalescent_tree.get_proband_descendants(unmapped_ancestor)
        variables = [self.variables[(proband, unmapped_ancestor)] for proband in unmapped_ancestor_proband_descendants]
        self.constraints.add(CommonAncestorConstraint(common_ancestor_coalescent_tree_id=unmapped_ancestor,
                                                      variables=variables))

    def add_intersect_constraints(self, unmapped_ancestor: int):
        unmapped_ancestor_probands = self.coalescent_tree.get_proband_descendants(unmapped_ancestor)
        for (a, b) in combinations(unmapped_ancestor_probands, 2):
            left_path_variable = self.variables[(a, unmapped_ancestor)]
            right_path_variable = self.variables[(b, unmapped_ancestor)]
            intersection_constraint = self.IntersectConstraint(left_path_variable, right_path_variable)
            self.constraints.add(intersection_constraint)

    def initialize_constraints(self):
        for level in self.coalescent_tree.levels[2:]:
            # Only the probands (vertices in the level 0) are mapped,
            # all the other vertices are to be mapped

            # The vertices on the first level cannot be used for building the "triplet constraints",
            # but we might use them for building other constraints
            for unmapped_ancestor in level:
                unmapped_ancestor_children: set = self.coalescent_tree.children_map[unmapped_ancestor]
                # Add the common ancestor constraint
                self.add_common_ancestor_constraint(unmapped_ancestor)
                # Add intersection constraints
                self.add_intersect_constraints(unmapped_ancestor)

    def get_pedigree_vertex_level(self, pedigree_vertex: int):
        return self.path_processed_graph.vertex_to_level_map[pedigree_vertex]

    def find_common_ancestor_related_variables(self, coalescent_vertex: id):
        return [self.variables[(proband, coalescent_vertex)] for proband
                in self.coalescent_tree.get_proband_descendants(coalescent_vertex)]

    def find_mapping(self):
        pass
