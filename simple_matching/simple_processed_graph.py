from genealogical_graph import GenealogicalGraph
from processed_graph import ProcessedGraph


class SimpleProcessedGraph(ProcessedGraph):

    def __init__(self, genealogical_graph: GenealogicalGraph):
        super().__init__(genealogical_graph)
        self.common_ancestry_map = dict()

    def initialize_vertex_common_ancestry_entry(self, vertex: int):
        self.common_ancestry_map[vertex] = dict([(x, {vertex}) for x in vertex.descendants])

    def update_vertex_common_ancestry_entry(self, vertex: int, common_ancestor_partner: int, parent: int):
        value = self.common_ancestry_map[parent][common_ancestor_partner]
        dictionary = self.common_ancestry_map[vertex]
        if common_ancestor_partner in self.common_ancestry_map[vertex]:
            dictionary[common_ancestor_partner].update(value)
        else:
            dictionary[common_ancestor_partner] = set(value)

    def get_common_ancestors_map(self, vertex: int):
        return self.common_ancestry_map[vertex]

    def align_with_coalescence_tree(self, coalescence_tree: GenealogicalGraph):
        sink_vertices_labels = [x for x in coalescence_tree.get_sink_vertices()]
        mapping = dict(zip(sink_vertices_labels, sink_vertices_labels))
        levels_to_map = coalescence_tree.levels[1:]
        for level in levels_to_map:
            for vertex in level:
                assert len(vertex.parents) < 2
                assert vertex.label not in mapping
                # Try to find the mapping for the vertex
                if len(vertex.children) < 2:
                    assert not vertex.parents
                    continue
                for child in vertex.children:
                    assert child in mapping
                ancestors: {int} = self.get_common_ancestors(vertex.children)
                ancestors = ancestors.difference(mapping.values())
                # Trying to assign every ancestor
                for ancestor in ancestors:
                    # Verifying that the ancestor is compatible with the coalescence tree
                    if not vertex.parents:
                        mapping[vertex.label] = ancestor
                        break
                    else:
                        [parent] = vertex.parents
                        parent = self.vertices[parent]
                        result = True
                        for descendant in parent.get_excluded_descendants(vertex):
                            if descendant in mapping:
                                common_ancestors = self.get_common_ancestors_pair(
                                    mapping[descendant], ancestor)
                                if not common_ancestors:
                                    print("Refusing mapping " + str(ancestor) + " for " + str(vertex.label))
                                    print()
                                    result = False
                                    break
                        if result:
                            mapping.update({vertex.label: ancestor})
                            break
                assert vertex.label in mapping
        return mapping

    def is_correct(self):
        vertices_number = len(self.vertices)
        indices = [x for x in self.vertices]
        indices.sort()
        assert indices == list(range(0, vertices_number))
        assert vertices_number == len(self.common_ancestry_map)
        for vertex, dictionary in self.common_ancestry_map.items():
            assert len(dictionary) <= vertices_number
            for key, value in dictionary.items():
                assert self.common_ancestry_map[key][vertex] == value

    def get_common_ancestors_pair(self, first, second):
        if second not in self.common_ancestry_map[first]:
            return []
        return self.common_ancestry_map[first][second]

    def get_common_ancestors(self, args: [int]):
        assert len(args) > 1
        [first, second] = args[:2]
        common_ancestors = self.get_common_ancestors_pair(first, second)
        for vertex in args[3:]:
            new_common_ancestors = set()
            for common_ancestor in common_ancestors:
                new_common_ancestors = \
                    new_common_ancestors.union(self.get_common_ancestors_pair(vertex, common_ancestor))
            common_ancestors = new_common_ancestors
        return common_ancestors

    def notify_processed_level(self, level_index: int):
        pass
