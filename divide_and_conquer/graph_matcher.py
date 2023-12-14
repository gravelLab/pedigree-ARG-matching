from divide_and_conquer.subtree_matcher import *


class GraphMather:

    def __init__(self, processed_graph: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree,
                 initial_mapping: dict = None):
        self.pedigree = processed_graph
        self.coalescent_tree = coalescent_tree
        if initial_mapping is None:
            initial_mapping = {x: x for x in processed_graph.probands}
        self.initial_mapping = initial_mapping

    def find_mapping(self):
        coalescent_tree_vertex_to_subtrees = dict()
        for vertex in self.coalescent_tree.levels[0]:
            coalescent_tree_vertex_to_subtrees[vertex] = get_subtree_matcher_for_coalescent_tree_proband(
                vertex, self.initial_mapping[vertex])
        for level in self.coalescent_tree.levels[1:]:
            for vertex in level:
                vertex_subtrees = get_subtrees_from_children(vertex, self.coalescent_tree, self.pedigree,
                                                             coalescent_tree_vertex_to_subtrees)
                coalescent_tree_vertex_to_subtrees[vertex] = vertex_subtrees
        # root_assignments = [coalescent_tree_vertex_to_subtrees[x] for x in self.coalescent_tree.levels[-1]]
        return coalescent_tree_vertex_to_subtrees
