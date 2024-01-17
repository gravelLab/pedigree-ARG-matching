"""!
@file descendant_memory_cache.py
@brief The file contains the basic implementation of the descendant cache that stores all the information in RAM.
"""


from descendant_cache import DescendantCache


class DescendantMemoryCache(DescendantCache):
    """
    The basic implementation of the descendant memory cache storing all the data in RAM.
    """

    def get_dictionary(self):
        return self.vertex_descendant_map

    def __init__(self):
        self.vertex_descendant_map = dict()

    def record_proband(self, proband: int):
        self.vertex_descendant_map[proband] = {proband}

    def record_child_descendants(self, parent_id: int, child_id: int):
        if child_id not in self.vertex_descendant_map:
            self.vertex_descendant_map[child_id] = {child_id}
        child_set: set = self.vertex_descendant_map[child_id]
        if parent_id in self.vertex_descendant_map:
            # All the next calls
            self.vertex_descendant_map[parent_id].update(child_set)
        else:
            # First time calling this function for parent_id
            self.vertex_descendant_map[parent_id] = child_set.union([parent_id])

    def get_vertex_descendants(self, vertex_id: int):
        return self.vertex_descendant_map.get(vertex_id)
