"""!
@file genealogical_graph.py
@brief This file contains the realization of the genealogical graph that assigns the vertices to their genealogical
levels and provides a level-by-level framework for preprocessing the graph and the utility class for processing
coalescent trees.
"""
from __future__ import annotations

from collections import defaultdict

import networkx as nx
from matplotlib import pyplot as plt

from graph.descendant_cache import DescendantCache
from graph.descendant_memory_cache import DescendantMemoryCache
from graph.simple_graph import *


class GenealogicalGraph(SimpleGraph):
    """!
    This class represents a genealogical graph that inherits from the /ref graph::Graph class and additionally
    assigns all the vertices within the graph to their levels. The level of a vertex is defined recursively as follows:
    1) All the proband individuals are assigned to level 0.
    2) If the maximal level among a vertex's children is n, then the vertex's level is n + 1.
    In other words, the level of a vertex is the length of the longest path from a proband to this vertex in a graph.
    """

    def __init__(self, pedigree: SimpleGraph = None, probands: [int] = None, initialize_levels: bool = True):
        if pedigree is None:
            super().__init__()
        else:
            self.parents_map = pedigree.parents_map
            self.children_map = pedigree.children_map
            self.vertices_number = pedigree.vertices_number
        if probands is None:
            self.initialize_probands()
        self.probands = probands
        self.levels = []
        self.vertex_to_level_map = dict()
        # TODO: Save the state of the levels (track whether or not they should be recalculated)
        if initialize_levels:
            self.initialize_vertex_to_level_map()
        self.descendant_writer = DescendantMemoryCache()

    def _find_cycle(self):
        visited = set()  # To track visited nodes
        recursion_stack = set()  # To track nodes currently in the recursion stack
        cycle = []  # To store the cycle if found

        def dfs(node, path):
            nonlocal cycle
            if node in recursion_stack:
                # A cycle is detected; capture the cycle path
                cycle_start_index = path.index(node)
                cycle = path[cycle_start_index:] + [node]
                return True

            if node in visited:
                return False

            visited.add(node)
            recursion_stack.add(node)
            path.append(node)

            for child in self.children_map.get(node, []):
                if dfs(child, path):
                    return True

            path.pop()
            recursion_stack.remove(node)
            return False

        for node in self.children_map:
            if node not in visited:
                if dfs(node, []):
                    break

        return cycle

    def initialize_vertex_to_level_map(self):
        """!
        @brief Assigns the vertices to the level. Refer to the docstring for this class to understand how a level of
        a vertex is defined.
        """
        self.initialize_probands()
        self.levels = []
        current_level = 1
        self.vertex_to_level_map = {x: 0 for x in self.probands}
        current_level_vertices = self.probands
        while current_level_vertices:
            current_level_vertices = {parent for x in current_level_vertices for parent in self.parents_map.get(x, [])}
            for vertex in current_level_vertices:
                self.vertex_to_level_map[vertex] = current_level
            current_level += 1
            self.levels.append(list())
        for vertex, level in self.vertex_to_level_map.items():
            self.levels[level].append(vertex)

    def initialize_probands(self):
        self.probands = {x for x in self.parents_map if not self.children_map.get(x, [])}

    def get_proband_descendants(self, vertex_id: int):
        """!
        @param vertex_id The vertex for which the proband descendants should be found.
        @return The descendants of the given vertex which are the probands.
        """
        return [x for x in self.get_vertex_descendants(vertex_id) if x in self.get_probands()]

    def set_descendant_writer(self, descendant_writer: DescendantCache):
        """!
        @brief The descendant writer to be set.
        """
        if self.descendant_writer is not None:
            raise Exception("Overriding the descendant writer is forbidden!")
        self.descendant_writer = descendant_writer

    @staticmethod
    def get_graph_from_tree(tree: Tree, probands=None):
        pedigree = SimpleGraph.get_graph_from_tree(tree)
        reformatted_parent_dict = defaultdict(list)
        for key, value in pedigree.parents_map.items():
            reformatted_parent_dict[key] = [value]
        pedigree.parents_map = reformatted_parent_dict
        graph = GenealogicalGraph(pedigree)
        return graph

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
        # Default behaviour goes here
        # Additional behaviour
        self.additional_process_level_vertex(parent_label, child, level)
        return parent_label

    def additional_process_level_vertex(self, parent_label: int, child: int, level: int):
        pass

    def get_vertex_descendants(self, vertex_id: int):
        return [x for x in self.descendant_writer.get_vertex_descendants(vertex_id) if x in self.vertex_to_level_map]

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

    def get_min_levels(self) -> (list, dict):
        min_levels_dict = dict()
        min_levels = []
        current_level = self.levels[0]
        level_index = 0
        while current_level:
            for vertex in current_level:
                min_levels_dict[vertex] = level_index
            min_levels.append(list(current_level))
            next_level = set()
            [next_level.update(self.parents_map[vertex]) for vertex in current_level if vertex in self.parents_map]
            current_level = next_level
            level_index += 1
        return min_levels, min_levels_dict

    def get_vertices_for_given_level(self, level):
        """!
        @brief Returns the vertices belonging to the specified level.
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
        if self.probands is None:
            self.initialize_probands()
        return self.probands

    def remove_edge(self, parent: int, child: int, recalculate_levels: bool = True):
        self.parents_map[child].remove(parent)
        self.children_map[parent].remove(child)
        if recalculate_levels:
            self.initialize_vertex_to_level_map()

    def add_edge(self, parent: int, child: int, recalculate_levels: bool = True):
        self.parents_map[child].append(parent)
        self.children_map[parent].append(child)
        if recalculate_levels:
            self.initialize_vertex_to_level_map()

    def remove_vertex(self, vertex: int, recalculate_levels: bool = True):
        """!
        @brief Removes the vertex from the graph.
        @param vertex The vertex to be removed.
        @param recalculate_levels: Specified whether the levels are to be recalculated after successfully removing
        the vertex from the graph
        """
        removed = super().remove_vertex(vertex)
        # Update the levels if they have been calculated
        if removed and self.probands:
            self.probands.discard(vertex)
            if recalculate_levels:
                self.initialize_vertex_to_level_map()
            else:
                vertex_level = self.vertex_to_level_map.pop(vertex)
                self.levels[vertex_level].remove(vertex)
        return removed

    def remove_vertices(self, vertices: [int], recalculate_levels: bool = True):
        """!
        @brief Removes the given vertices from the graph.
        @param vertices The vertices to be removed.
        @param recalculate_levels: Specified whether the levels are to be recalculated after removing the vertices
        from the graph
        """
        # vertices_to_remove = frozenset(vertices)
        # for vertex in vertices_to_remove:
        #     self.children_map.pop(vertex, None)
        #     self.parents_map.pop(vertex, None)
        # for child, current_parents in list(self.parents_map.items()):
        #     self.parents_map[child] = [x for x in current_parents if x not in vertices_to_remove]
        # for parent, current_children in list(self.children_map.items()):
        #     self.children_map[parent] = [x for x in current_children if x not in vertices_to_remove]
        # if recalculate_levels:
        #     self.initialize_vertex_to_level_map()
        # else:
        #     for vertex in vertices_to_remove:
        #         self.vertex_to_level_map.pop(vertex, None)
        #     self.levels = [
        #         [x for x in level if x not in vertices_to_remove] for level in self.levels
        #     ]
        for vertex in vertices:
            self.remove_vertex(vertex, False)
        if recalculate_levels:
            self.initialize_vertex_to_level_map()

    def get_ascending_genealogy_from_vertices_by_levels(self, vertices: [int]):
        ascending_genealogy = self.get_ascending_genealogy_from_vertices(vertices)
        level_to_ascending_vertices = defaultdict(list)
        # Group vertices in ascending_genealogy based on their levels
        for vertex in ascending_genealogy:
            level = self.vertex_to_level_map[vertex]
            level_to_ascending_vertices[level].append(vertex)
        ascending_genealogy_by_levels = [
            level_to_ascending_vertices.get(i, []) for i in range(len(self.levels))
        ]
        assert len(ascending_genealogy) == sum(len(sublist) for sublist in ascending_genealogy_by_levels)
        return ascending_genealogy_by_levels

    def write_levels_to_file(self, file, levels):
        processed_individuals = set()
        for level in levels:
            for vertex_ploid_id in level:
                vertex_individual_id = vertex_ploid_id // 2
                if vertex_individual_id in processed_individuals:
                    continue
                processed_individuals.add(vertex_individual_id)
                [first_parent_id, second_parent_id] = [-1, -1]
                if vertex_ploid_id in self.parents_map:
                    ploid_id = 2 * vertex_individual_id
                    ploid_parents = self.parents_map.get(ploid_id, [])
                    if ploid_parents:
                        [first_parent, _] = ploid_parents
                        first_parent_id = first_parent // 2
                    ploid_id += 1
                    ploid_parents = self.parents_map.get(ploid_id, [])
                    if ploid_parents:
                        [second_parent, _] = ploid_parents
                        second_parent_id = second_parent // 2
                file.write(f"{vertex_individual_id} {first_parent_id} {second_parent_id}\n")

    def draw_graph(self):
        g = nx.DiGraph()
        for level_index, level in enumerate(self.levels):
            for vertex in level:
                g.add_node(vertex, time=level_index)
                vertex_parents = self.parents_map.get(vertex, [])
                for parent in vertex_parents:
                    g.add_edge(vertex, parent)
        fig, ax = plt.subplots(figsize=(655, 655))
        node_size = 1000
        pos = nx.multipartite_layout(g, subset_key="time", align="horizontal")
        nx.draw_networkx(g, pos, ax=ax, node_size=node_size, with_labels=True)
        plt.savefig("small_graph.png")

    def save_ascending_genealogy_to_file(self, filepath: str, probands: [int]):
        levels = self.get_ascending_genealogy_from_vertices_by_levels(probands)
        file = open(filepath, 'w')
        self.write_levels_to_file(file, levels)
        file.close()

    def save_to_file(self, filepath: str):
        file = open(filepath, 'w')
        self.write_levels_to_file(file, self.levels)
        file.close()

    def reduce_to_ascending_genealogy(self, probands: [int], recalculate_levels: bool = True):
        ascending_genealogy = self.get_ascending_genealogy_from_vertices(probands)
        graph_vertices = self.get_vertices()
        vertices_to_remove = graph_vertices.difference(ascending_genealogy)
        self.remove_vertices(vertices=vertices_to_remove, recalculate_levels=recalculate_levels)

    @staticmethod
    def get_diploid_graph_from_file(filepath: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False, initialize_levels: bool = True) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a diploid genealogical graph from a file.
        @param filepath The filename from which the graph will be read.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        Refer to the documentation for SimpleGraph.get_graph_from_file
        @param separation_symbol The separation sequence used in the file. Refer to the documentation for
        SimpleGraph.get_graph_from_file
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @returns The parsed graph.
        @param initialize_levels: Specifies whether the levels should be initialized after parsing
        """
        pedigree: SimpleGraph = SimpleGraph.get_diploid_graph_from_file(filepath=filepath,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree, initialize_levels=initialize_levels)

    @staticmethod
    def get_haploid_graph_from_file(filepath: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a haploid genealogical graph from a file.
        @param filepath The filename from which the graph will be read.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        Refer to the documentation for SimpleGraph.get_graph_from_file
        @param separation_symbol The separation sequence used in the file. Refer to the documentation for
        SimpleGraph.get_graph_from_file
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @returns The parsed graph.
        """
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filepath,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree)
