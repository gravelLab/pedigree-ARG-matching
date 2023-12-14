from descendant_cache import DescendantCache
from genealogical_graph import GenealogicalGraph


class DescendantNoMemoryCache(DescendantCache):

    def __init__(self, genealogical_graph: GenealogicalGraph):
        self.genealogical_graph = genealogical_graph

    def record_proband(self, proband: int):
        pass

    def get_vertex_descendants(self, vertex_id: int):
        result = {vertex_id}
        for child in self.genealogical_graph.children_map[vertex_id]:
            result.update(self.get_vertex_descendants(child))
        return result

    def get_dictionary(self):
        return dict()

    def record_child_descendants(self, parent_id, child_id):
        pass
