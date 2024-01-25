"""!
@file genealogical_graph.py
@brief This file contains the realization of the genealogical graph that assigns the vertices to their genealogical
levels and provides a level-by-level framework for preprocessing the graph and the utility class for processing
coalescent trees.
"""

import itertools

from descendant_cache import DescendantCache
from descendant_memory_cache import DescendantMemoryCache
from graph import *


class GenealogicalGraph(Graph):
    """!
    This class represents a genealogical graph that inherits from the /ref graph::Graph class and additionally
    assigns all the vertices within the graph to their levels. The level of a vertex is defined recursively as follows:
    1) All the proband individuals are assigned to level 0.
    2) If the maximal level among a vertex's children is n, then the vertex's level is n + 1.
    In other words, the level of a vertex is the length of the longest path from a proband to this vertex in a graph.
    """

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
        """!
        @brief Assigns the vertices to the level. Refer to the docstring for this class to understand how a level of
        a vertex is defined.
        """
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
        """!
        @param vertex_id The vertex for which the proband descendants should be found.
        @return The descendants of the given vertex which are the probands.
        """
        # TODO: Cache the list of probands
        return [x for x in self.get_vertex_descendants(vertex_id) if x in self.get_probands()]

    def set_descendant_writer(self, descendant_writer: DescendantCache):
        """!
        @brief The descendant writer to be set.
        """
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
        """!
        @bried Returns the vertices belonging to the specified level.
        @param level The level to be used.
        """
        return self.levels[level]

    def get_top_level_vertices(self):
        """!
        @brief Returns the vertices at the top level of the graph.
        """
        return self.levels[-1]

    def get_probands(self):
        """!
        @brief Returns the graph's probands.
        """
        return self.probands

    def remove_vertex(self, vertex: int):
        """!
        @brief Removes the vertex from the graph.
        @param vertex The vertex to be removed.
        """
        super().remove_vertex(vertex)
        # Vertex usually must belong to this map, this can only happen if we try to remove the same vertex twice
        if vertex in self.vertex_to_level_map:
            vertex_level = self.vertex_to_level_map[vertex]
            self.vertex_to_level_map.pop(vertex)
            level: [int] = self.levels[vertex_level]
            level.remove(vertex)


class CoalescentTree(GenealogicalGraph):
    """!
    This is a helper class that is responsible for working with coalescent trees. Apart from the functionality
    of the GenealogicalGraph, it calculates the connected components (clades) of the graph.
    """

    def __init__(self, pedigree: Graph):
        super().__init__(pedigree=pedigree)
        self.descendant_writer = DescendantMemoryCache()
        self.initialize_genealogical_graph_from_probands()
        self.find_clades()

    def find_clades(self):
        """!
        TODO: Add a method that stores the clades in the coalescent tree after adding the vertices set to Graph.
        """
        pass

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


    # @staticmethod
    # def get_coalescent_tree(tree: Tree, probands=None):
    #     genealogical_graph = GenealogicalGraph.get_graph_from_tree(tree, probands)
    #     return CoalescentTree(genealogical_graph)
