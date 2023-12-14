from enum import Enum
from io import TextIOWrapper
from genealogical_graph import GenealogicalGraph
from descendant_cache import DescendantCache

default_file_prefix = ""
default_file_postfix = ".dat"


def get_file_to_write_for_vertex(vertex_id: int):
    return open(f"{default_file_prefix}{vertex_id}{default_file_postfix}", "a")


def read_level_from_disk(level: int):
    filename = f"{default_file_prefix}{level}{default_file_postfix}"
    # TODO Open file and read the descendants
    return []


class CacheState(Enum):
    READ = 0,
    WRITE = 1


class CacheMemoryMode(Enum):
    FRUGAL = 0,
    FAST = 1


class DescendantFileCache(DescendantCache):

    def __init__(self, graph: GenealogicalGraph, memory_mode: CacheMemoryMode = CacheMemoryMode.FAST):
        self.graph = graph
        graph.set_descendant_writer(self)

        # Map the vertex id to the line at which it is stored
        level_indices = range(0, len(graph.levels))
        self.vertex_line_map = dict()
        self.vertex_file_map = {x: get_file_to_write_for_vertex(x) for x in graph.}
        self.file_current_line_map = {x: 0 for x in level_indices}

        self.cache_high_level = dict()
        self.cache_low_level = {x: [x] for x in graph.get_probands()}
        self.current_level = 1
        self.line_count = 0
        self.current_state = CacheState.READ
        self.memory_mode = memory_mode

    def switch_to_write(self):
        self.current_state = CacheState.WRITE
        for file in self.file_map:
            file: TextIOWrapper
            file.close()
        del self.file_current_line_map

    def write_to_disk(self, vertex_id: int, descendant_list: [int]):
        # Keep a list of all the files in the class and close them in the switch_to_write method
        # Open the corresponding file and append the descendant list to the end of the file
        # Record the line which keeps the descendant list of this vertex
        level: int = self.graph.vertex_to_level_map[vertex_id]
        file: TextIOWrapper = self.file_map.get(level)
        file.write(descendant_list)
        file.write('\n')
        current_line: int = self.file_current_line_map.get(level)
        self.file_current_line_map[level] = current_line + 1
        self.vertex_line_map[vertex_id] = current_line

    # Assuming that child belongs to the low level and parents belong to the high level,
    # update the left parent and right parent's entries with child's descendants
    # Write child's descendants to the disk
    def record_descendants(self, left_parent_id: int, right_parent_id: int, child_id: int):
        assert self.current_state == CacheState.READ
        self.vertex_line_map.update({left_parent_id: self.line_count})
        self.line_count += 1
        self.vertex_line_map.update({right_parent_id: self.line_count})
        child_descendants = self.cache_low_level[child_id]
        if self.memory_mode == CacheMemoryMode.FRUGAL:
            self.cache_low_level.pop(child_id)
        if left_parent_id not in self.cache_high_level:
            self.cache_high_level[left_parent_id] = list(child_descendants)
        else:
            self.cache_high_level[left_parent_id].extend()

    def get_vertex_descendants(self, vertex_id: int):
        assert self.current_level == CacheState.WRITE
        pass

    def climb_level(self):
        # Update the index
        self.current_level += 1
        assert self.current_level < len(self.graph.levels)
        # Reassign the lower level
        self.cache_low_level = self.cache_high_level
        if self.current_state == CacheState.READ:
            # Close the file and open a new one
            self.current_file.close()
            self.current_file = get_file_to_write_for_level(self.current_level)
            self.line_count = 0
        else:
            self.cache_high_level = read_level_from_disk(self.current_level)

    def descend_level(self):
        self.current_level -= 1
        assert self.current_level > -1
        self.cache_high_level = self.cache_low_level
        self.cache_low_level = read_level_from_disk(self.current_level - 1)
