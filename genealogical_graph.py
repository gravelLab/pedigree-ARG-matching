"""!
@file genealogical_graph.py
@brief This file contains the realization of the genealogical graph that assigns the vertices to their genealogical
levels and provides a level-by-level framework for preprocessing the graph and the utility class for processing
coalescent trees.
"""
from __future__ import annotations
from tskit import TreeSequence

import networkx as nx
from matplotlib import pyplot as plt

from descendant_cache import DescendantCache
from descendant_memory_cache import DescendantMemoryCache
from simple_graph import *


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
            probands = {x for x in self.parents_map if x not in self.children_map}
        self.probands = probands
        self.levels = []
        self.vertex_to_level_map = dict()
        if initialize_levels:
            self.initialize_vertex_to_level_map()
        self.descendant_writer = DescendantMemoryCache()

    def initialize_vertex_to_level_map(self):
        """!
        @brief Assigns the vertices to the level. Refer to the docstring for this class to understand how a level of
        a vertex is defined.
        """
        self.levels = []
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

    def remove_vertex(self, vertex: int, recalculate_levels: bool = True):
        """!
        @brief Removes the vertex from the graph.
        @param vertex The vertex to be removed.
        @param recalculate_levels: Specified whether the levels are to be recalculated after successfully removing
        the vertex from the graph
        """
        removed = super().remove_vertex(vertex)
        # Vertex usually must belong to this map, this can only happen if we try to remove the same vertex twice
        if removed:
            vertex_level = self.vertex_to_level_map[vertex]
            self.vertex_to_level_map.pop(vertex)
            level: [int] = self.levels[vertex_level]
            level.remove(vertex)
            if vertex in self.probands:
                self.probands.remove(vertex)
            if recalculate_levels:
                self.initialize_vertex_to_level_map()
        return removed

    def remove_vertices(self, vertices: [int], recalculate_levels: bool = True):
        """!
        @brief Removes the given vertices from the graph.
        @param vertices The vertices to be removed.
        @param recalculate_levels: Specified whether the levels are to be recalculated after removing the vertices
        from the graph
        """
        for vertex in vertices:
            self.remove_vertex(vertex, False)
        if recalculate_levels:
            self.initialize_vertex_to_level_map()

    def get_ascending_genealogy_from_vertices_by_levels(self, vertices: [int]):
        ascending_genealogy = self.get_ascending_genealogy_from_vertices(vertices)
        return [[x for x in level if x in ascending_genealogy] for level in self.levels]

    def write_levels_to_file(self, file, levels):
        processed_ids = set()
        for level in levels:
            for vertex in level:
                vertex_id = vertex // 2
                if vertex_id in processed_ids:
                    continue
                else:
                    processed_ids.add(vertex_id)
                [first_parent_id, second_parent_id] = [-1, -1]
                if vertex in self.parents_map:
                    ploid_id = 2 * vertex_id
                    if ploid_id in self.parents_map:
                        [first_parent, _] = self.parents_map[ploid_id]
                        first_parent_id = first_parent // 2
                    ploid_id += 1
                    if ploid_id in self.parents_map:
                        [second_parent, _] = self.parents_map[ploid_id]
                        second_parent_id = second_parent // 2
                file.write(f"{vertex_id} {first_parent_id} {second_parent_id}\n")

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

    def save_ascending_genealogy_to_file(self, filename: str, probands: [int]):
        levels = self.get_ascending_genealogy_from_vertices_by_levels(probands)
        file = open(filename, 'w')
        self.write_levels_to_file(file, levels)
        file.close()

    def save_to_file(self, filename: str):
        file = open(filename, 'w')
        self.write_levels_to_file(file, self.levels)
        file.close()

    @staticmethod
    def get_diploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a diploid genealogical graph from a file.
        @param filename The filename from which the graph will be read.
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
        pedigree: SimpleGraph = SimpleGraph.get_diploid_graph_from_file(filename=filename,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree)

    @staticmethod
    def get_haploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a haploid genealogical graph from a file.
        @param filename The filename from which the graph will be read.
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
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filename,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree)


class CoalescentTree(GenealogicalGraph):
    """!
    This is a helper class that is responsible for working with coalescent trees. Apart from the functionality
    of the GenealogicalGraph, it calculates the connected components (clades) of the graph.
    """

    def __init__(self, pedigree: SimpleGraph):
        super().__init__(pedigree=pedigree)
        self.descendant_writer = DescendantMemoryCache()
        self.initialize_genealogical_graph_from_probands()

    def get_largest_clade_by_size(self):
        clades = self.get_connected_components()
        largest_clade = max(clades, key=len)
        return largest_clade

    def get_largest_clade_by_probands(self):
        def intersection_size(clade):
            return len(self.probands.intersection(clade))

        clades = self.get_connected_components()
        largest_clade = max(clades, key=lambda clade: intersection_size(clade))
        return largest_clade

    def get_root_for_clade(self, clade: [int]):
        max_level_vertex = max(clade, key=lambda x: self.vertex_to_level_map[x])
        max_level = self.vertex_to_level_map[max_level_vertex]
        root_vertices = [x for x in clade if x in self.levels[max_level]]
        if len(root_vertices) != 1:
            raise Exception("Invalid clade value")
        assert root_vertices[0] == max_level_vertex
        return root_vertices[0]

    def get_roots_for_clade(self, clade: [int]):
        max_level_vertex = max(clade, key=lambda x: self.vertex_to_level_map[x])
        max_level = self.vertex_to_level_map[max_level_vertex]
        return [x for x in clade if x in self.levels[max_level]]

    def get_subtree_from_vertices(self, vertices: [int]):
        def narrow_function(map_to_narrow: dict):
            return {key: map_to_narrow[key] for key in vertices if key in map_to_narrow}

        children_map_vertices = narrow_function(self.children_map)
        parents_map_vertices = narrow_function(self.parents_map)
        vertices_number = len(vertices)
        pedigree = SimpleGraph(children_map=defaultdict(list, children_map_vertices),
                               parents_map=defaultdict(list, parents_map_vertices),
                               vertices_number=vertices_number)
        return CoalescentTree(pedigree=pedigree)

    def remove_unary_nodes(self):
        """
        Removes all the unary nodes in the coalescent tree and recalculates the levels of the coalescent tree.
        """
        for level in self.levels[1:].__reversed__():
            intermediate_nodes = []
            for vertex in level:
                # Since the first level is omitted, all the vertices processed here must have children
                children = self.children_map[vertex]
                if len(children) == 1:
                    child = vertex
                    while len(children) == 1:
                        intermediate_nodes.append(child)
                        [child] = children
                        if child not in self.children_map:
                            break
                        children = self.children_map[child]
                    if vertex in self.parents_map:
                        [parent] = self.parents_map[vertex]
                        self.children_map[parent].append(child)
                        self.parents_map[child] = [parent]
            self.remove_vertices(intermediate_nodes, False)
        self.initialize_vertex_to_level_map()
        assert not [x for x in self.children_map if len(self.children_map[x]) == 1]

    # def save_to_file(self, filename: str):
    #     file = open(filename, 'w')
    #     self.write_levels_to_file(file, self.levels)
    #     file.close()
    #
    # def save_ascending_genealogy_to_file(self, filename: str, probands: [int]):
    #     levels = self.get_ascending_genealogy_from_vertices_by_levels(probands)
    #     file = open(filename, 'w')
    #     self.write_levels_to_file(filename, levels)
    #     file.close()

    def write_levels_to_file(self, file, levels):
        for level in levels:
            for vertex in level:
                if vertex in self.parents_map:
                    parents = self.parents_map[vertex]
                    file.write(f"{vertex} {' '.join(str(parent) for parent in parents)}\n")

    @staticmethod
    def get_coalescent_tree_from_file(filename: str, max_parent_number: int = 2,
                                      missing_parent_notation=None, separation_symbol=' ',
                                      skip_first_line: bool = False):
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filename,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        result = CoalescentTree(pedigree=pedigree)
        return result

    @staticmethod
    def get_coalescent_tree(tree: Tree, probands=None):
        genealogical_graph = GenealogicalGraph.get_graph_from_tree(tree, probands)
        return CoalescentTree(genealogical_graph)

    @staticmethod
    def get_arg(tree_sequence: TreeSequence):
        first_tree = tree_sequence.first()
        coalescent_tree = CoalescentTree.get_coalescent_tree(first_tree)
        present_edges = {(key, value) for key, value in first_tree.parent_dict.items()}
        for tree in tree_sequence.trees():
            for (child, parent) in tree.parent_dict.items():
                if not (child, parent) in present_edges:
                    coalescent_tree.children_map[parent].append(child)
                    coalescent_tree.parents_map[child].append(parent)
        return coalescent_tree
