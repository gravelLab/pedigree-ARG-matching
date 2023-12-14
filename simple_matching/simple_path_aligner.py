from path_matching.path_processed_graph import PathProcessedGraph, PathAwareGenealogicalGraph
from genealogical_graph import GenealogicalGraph, CoalescentTree
from functools import cmp_to_key


# Create an iterator that will return the vertices one by one. It should have a constructor accepting a genealogical
# graph and a copy method that will simply copy a pointer and an integer

# Use this iterator to perform the recursive calls. Whenever you fail, you will have the previous assignment state

class AlignmentIterator:

    def __init__(self, vertices: [int], position: int = 0):
        self.vertices = vertices
        self.position = position
        self.vertices_length = len(vertices)

    @staticmethod
    def get_from_genealogical_graph(genealogical_graph: GenealogicalGraph):
        # Since currently a vertex can be present in different levels, we should find the maximum level for every vertex
        vertex_to_level = dict()
        for index in range(1, len(genealogical_graph.levels)):
            for vertex in genealogical_graph.levels[index]:
                vertex_to_level[vertex.label] = index

        def compare_vertices(first: int, second: int):
            return vertex_to_level[first] - vertex_to_level[second]

        result = sorted(vertex_to_level.keys(), key=cmp_to_key(compare_vertices))
        return AlignmentIterator(vertices=result)

    def copy(self):
        return AlignmentIterator(self.vertices, self.position)

    def next(self):
        result = self.vertices[self.position]
        self.position = self.position + 1
        return result

    def finished(self):
        return self.vertices_length == self.position + 1


class SimplePathAligner:

    def __init__(self, path_graph: PathProcessedGraph, coalescence_tree: CoalescentTree,
                 parent_alignment=None, partial_mapping: dict = None, iterator: AlignmentIterator = None):
        self.path_graph = path_graph
        self.coalescence_tree = coalescence_tree
        self.parent_alignment = parent_alignment
        if partial_mapping is None:
            partial_mapping = dict()
            for proband in coalescence_tree.levels[0]:
                partial_mapping[proband.label] = proband.label
        self.partial_mapping = partial_mapping
        if iterator is None:
            iterator = AlignmentIterator.get_from_genealogical_graph(self.path_graph)
        self.iterator = iterator

    def align(self):
        iterator = self.iterator
        while not iterator.finished():
            # TODO Initialize the partial mapping with None values
            # Reading the next vertex to be mapped
            vertex_to_map = iterator.next()
            # Getting the candidate vertices for the assignment
            candidates = self.get_candidates_for_vertex(vertex_to_map)
            if candidates is None:
                self.partial_mapping[vertex_to_map] = None
                assert vertex_to_map not in self.coalescence_tree.parents_map
                continue
            result = list()
            for candidate in candidates:
                # TODO Redesign this function, so that it returns a list of PathAligner objects
                if self.verify_assignment(vertex_to_map=vertex_to_map, candidate=candidate):
                    result.append(candidate)
            if not result:
                return None
            while len(result) > 1:
                candidate = result.pop()
                copy_mapping = dict(self.partial_mapping)
                copy_mapping[vertex_to_map] = candidate
                copy_iterator = iterator.copy()
                aligner = SimplePathAligner(path_graph=self.path_graph, coalescence_tree=self.coalescence_tree,
                                            parent_alignment=self, partial_mapping=copy_mapping, iterator=copy_iterator)
                alignment_result = aligner.align()
                if alignment_result is not None:
                    return alignment_result
            self.partial_mapping[vertex_to_map] = result.pop()
        return self.partial_mapping

    def get_candidates_for_vertex(self, coalescence_tree_vertex: int):
        if coalescence_tree_vertex not in self.coalescence_tree.children_map:
            return None
        coalescence_tree_children = self.coalescence_tree.children_map[coalescence_tree_vertex]
        for child in coalescence_tree_children:
            assert child in self.partial_mapping
        if len(coalescence_tree_children) < 2:
            return None
        result = self.path_graph.get_common_ancestors(
            [self.partial_mapping[x] for x in coalescence_tree_children])
        return {key for key in result if key not in self.partial_mapping.values()}

    def verify_assignment(self, vertex_to_map: int, candidate: int):
        if vertex_to_map not in self.coalescence_tree.parents_map:
            return True
        assert vertex_to_map in self.coalescence_tree.children_map

        vertex_to_map_parent = self.coalescence_tree.parents_map[vertex_to_map]
        outgroup_descendants = [self.partial_mapping[x] for x in self.coalescence_tree.get_excluded_proband_descendants(
            parent_vertex=vertex_to_map_parent, child_vertex=vertex_to_map
        )]
        if not outgroup_descendants:
            return candidate in self.path_graph.parents_map
        vertex_to_map_descendants = [self.partial_mapping[x] for x in
                                     self.coalescence_tree.proband_descendants[vertex_to_map]]
        for descendant in vertex_to_map_descendants:
            paths_from_child_to_candidate = (self.path_graph.path_common_ancestry_map[descendant][candidate][candidate]
                                             .get_paths_for_vertex(descendant))
            for outgroup_vertex in outgroup_descendants:
                for common_ancestry in self.path_graph.path_common_ancestry_map[outgroup_vertex][descendant].values():
                    common_ancestry: PathProcessedGraph.PathCommonAncestry
                    descendant_coalescence_grandparent_paths = common_ancestry.get_paths_for_vertex(descendant)
                    # Verify that the assignment is consistent
                    for grandparent_path in descendant_coalescence_grandparent_paths:
                        for path_to_candidate in paths_from_child_to_candidate:
                            grandparent_path: PathAwareGenealogicalGraph.Path
                            path_to_candidate: PathAwareGenealogicalGraph.Path
                            if not grandparent_path.path.endswith(path_to_candidate.path):
                                return True
        return False
