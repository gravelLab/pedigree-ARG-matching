import itertools
import sys

from descendant_cache import DescendantCache
from descendant_memory_cache import DescendantMemoryCache
from graph import *


class GenealogicalGraph(Graph):
    def __init__(self, pedigree: Graph = None, probands: [int] = None, initialize_levels: bool = True):
        self.descendant_writer = None
        if pedigree is None:
            super().__init__()
        else:
            self.parents_map = pedigree.parents_map
            self.children_map = pedigree.children_map
            self.sink_vertices = pedigree.sink_vertices
            self.vertices_number = pedigree.vertices_number
        self.probands = probands
        self.levels = []
        self.vertex_to_level_map = dict()
        if initialize_levels:
            self.initialize_vertex_to_level_map()

    def initialize_vertex_to_level_map(self):
        if self.probands is None:
            self.probands = {x for x in self.parents_map if x not in self.children_map}
        current_level = 1
        self.vertex_to_level_map = {x: 0 for x in self.probands}
        current_level_vertices = self.probands
        while current_level_vertices:
            current_level_vertices = set(itertools.chain.from_iterable(
                [self.parents_map[x] for x in current_level_vertices if x in self.parents_map])
            )
            for vertex in current_level_vertices:
                self.vertex_to_level_map[vertex] = current_level
            current_level += 1
            self.levels.append(list())
        for vertex, level in self.vertex_to_level_map.items():
            self.levels[level].append(vertex)

    def get_proband_descendants(self, vertex_id: int):
        # TODO: Cache the list of probands
        return [x for x in self.get_vertex_descendants(vertex_id) if x in self.get_probands()]

    def set_descendant_writer(self, descendant_writer: DescendantCache):
        if self.descendant_writer is not None:
            raise Exception("Overriding the descendant writer is forbidden!")
        self.descendant_writer = descendant_writer

    # @staticmethod
    # def get_graph_from_tree(tree: Tree, probands=None):
    #     pedigree = Graph.get_graph_from_tree(tree)
    #     reformatted_parent_dict = dict()
    #     for key, value in pedigree.parents_map.items():
    #         reformatted_parent_dict[key] = [value]
    #     pedigree.parents_map = reformatted_parent_dict
    #     graph = GenealogicalGraph(pedigree)
    #     return graph

    def process_proband_vertex(self, proband_label: int):
        # Default behaviour
        self.descendant_writer.record_proband(proband_label)
        # Additional processing
        self.additional_process_proband_vertex(proband_label)

    def get_vertex_excluded_descendants(self, vertex_id: int, excluded_vertex: int):
        result = set()
        for child in self.children_map[vertex_id]:
            if child != excluded_vertex:
                result.update(self.get_vertex_descendants(child))
        return result

    def additional_process_proband_vertex(self, proband_label: int):
        pass

    def process_level_vertex(self, parent_label: int, child: int, level: int):
        # Default behaviour
        pass
        # Additional behaviour
        self.additional_process_level_vertex(parent_label, child, level)
        return parent_label

    def additional_process_level_vertex(self, parent_label: int, child: int, level: int):
        pass

    def get_vertex_descendants(self, vertex_id: int):
        return self.descendant_writer.get_vertex_descendants(vertex_id)

    def initialize_genealogical_graph_from_probands(self):
        for x in self.probands:
            self.process_proband_vertex(x)
        counter = 1
        for current_level_vertices in self.levels[1:]:
            for parent in current_level_vertices:
                if parent in self.children_map:
                    children = self.children_map[parent]
                    for child in children:
                        self.process_level_vertex(parent, child, counter)
                        self.descendant_writer.record_child_descendants(parent_id=parent, child_id=child)
            counter += 1
        return self

    def get_vertices_for_given_level(self, level):
        return self.levels[level]

    def get_top_level_vertices(self):
        return self.levels[-1]

    def get_probands(self):
        return self.probands

    def calculate_vertices_in_descendant_writer(self):
        integersMemory = 0
        setMemory = 0
        count = 0
        for vertex in self.vertex_to_level_map.keys():
            descendants = self.descendant_writer.get_vertex_descendants(vertex)
            if descendants is not None:
                count += len(descendants)
                setMemory += sys.getsizeof(descendants)
                for descendant in descendants:
                    integersMemory += sys.getsizeof(descendant)
        print(f"Integers occupy {integersMemory / (1024 * 1024)} MB of memory")
        print(f"Average integer take {integersMemory / count}")
        print(f"Total set memory: {setMemory}")
        return count

    def remove_vertex(self, vertex: int):
        super().remove_vertex(vertex)
        # Vertex usually must belong to this map, this can only happen if we try to remove the same vertex twice
        if vertex in self.vertex_to_level_map:
            vertex_level = self.vertex_to_level_map[vertex]
            self.vertex_to_level_map.pop(vertex)
            level: [int] = self.levels[vertex_level]
            level.remove(vertex)


class CoalescentTree(GenealogicalGraph):

    def __init__(self, pedigree: Graph):
        super().__init__(pedigree=pedigree)
        self.descendant_writer = DescendantMemoryCache()
        self.initialize_genealogical_graph_from_probands()
        self.proband_descendants = dict()
        for vertex in self.children_map.keys():
            vertex_proband_descendants = [x for x in self.get_vertex_descendants(vertex) if x in self.probands]
            self.proband_descendants.update({vertex: vertex_proband_descendants})
        self.verify_correctness()

    def remove_non_coalescing_nodes(self):
        modified_children_map = dict(self.children_map)
        for vertex, children in modified_children_map.items():
            vertex: int
            if len(children) == 0:
                raise Exception("Invalid children map")
            if len(children) == 1:
                # Remove the whole subgraph
                vertices = self.get_connected_component_for_vertex(vertex)
                for connected_component_vertex in vertices:
                    self.remove_vertex(connected_component_vertex)

    def get_excluded_proband_descendants(self, parent_vertex: int, child_vertex: int):
        children = list(self.children_map[parent_vertex])
        children: [int]
        children.remove(child_vertex)
        for x in children:
            if x not in self.proband_descendants:
                raise Exception(F"Vertex {x} is not in the proband descendants map!")
        return list(itertools.chain.from_iterable([self.proband_descendants[x] for x in children]))

    @staticmethod
    def get_coalescent_tree_from_file(filename):
        pedigree: Graph = Graph()
        file = open(filename)
        for line in file.readlines():
            (child, parent) = list(map(lambda name: int(name), line.strip('\n').split(' ')))
            pedigree.add_child(parent=parent, child=child)
            pedigree.parents_map[child] = [parent]
        file.close()
        result = CoalescentTree(pedigree=pedigree)
        return result

    def verify_correctness(self):
        pass
        # for level in self.levels[1:]:
        #     for vertex in level:
        #         assert len(self.get_proband_descendants(vertex)) > 1

    # @staticmethod
    # def get_coalescent_tree(tree: Tree, probands=None):
    #     genealogical_graph = GenealogicalGraph.get_graph_from_tree(tree, probands)
    #     return CoalescentTree(genealogical_graph)
