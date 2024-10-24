"""!
@file graph.py
@brief The file contains the realization of the basic graph class that is responsible for parsing a file with a
genealogical graph (i.e. a pedigree or a coalescent tree) and building the parents and children maps for this graph.
"""
from __future__ import annotations
import itertools
import warnings

from tskit import Tree
from collections import defaultdict


class SimpleGraph:
    """!
    This class represents a simple genealogical graph and stores the information about the vertices and the edges
    (children-parent relationships) in a graph
    """

    def __init__(self, children_map: dict = None, parents_map: dict = None, vertices_number: int = 0):
        if parents_map is None:
            parents_map = dict()
        if children_map is None:
            children_map = dict()
        self.children_map = children_map
        self.parents_map = parents_map
        self.vertices_number = vertices_number

    def get_connected_components(self):
        """!
        This method returns the connected components of the graph.
        """
        visited = set()
        connected_components = []
        for vertex in self.parents_map:
            if vertex in visited:
                continue
            connected_component = self.get_connected_component_for_vertex(vertex)
            assert not [x for x in connected_component if x in visited]
            visited.update(connected_component)
            connected_components.append(connected_component)
        return connected_components

    def get_connected_component_for_vertex(self, vertex: int):
        """!
        @brief The function finds all the vertices that are in the same connected component as the passed vertex.
        @param vertex The vertex for which the connected component should be calculated.
        @return The set of vertices that are in the same connected component as the passed vertex.
        """
        visited = set()
        component = []

        def dfs(source_vertex):
            if source_vertex not in visited:
                visited.add(source_vertex)
                component.append(source_vertex)
                for neighbour in self.parents_map.get(source_vertex, []):
                    dfs(neighbour)
                for neighbour in self.children_map.get(source_vertex, []):
                    dfs(neighbour)

        dfs(vertex)
        return component

    def get_minimal_path_length(self, descendant: int, ancestor: int) -> int:
        current_level_vertices = {descendant}
        length = 0
        while ancestor not in current_level_vertices:
            length += 1
            current_level_vertices = set(itertools.chain.from_iterable(
                [self.parents_map[x] for x in current_level_vertices if x in self.parents_map])
            )
        return length

    def get_ascending_genealogy_from_vertices(self, vertices: [int]):
        """!
        @brief This method returns all the vertices in the ascending genealogy for the given list of vertices.
        @param vertices The vertices for which the ascending genealogy should be calculated.
        """
        result = set(vertices)
        current_level_vertices = set(result)
        while current_level_vertices:
            next_level_vertices = set()
            for vertex in current_level_vertices:
                parents = self.parents_map.get(vertex, [])
                result.update(parents)
                next_level_vertices.update(parents)
            current_level_vertices = next_level_vertices
        return result

    def remove_vertex(self, vertex: int):
        """!
        @brief This function removes the vertex from the graph.
        @param vertex The vertex to be removed.
        @return Returns whether the vertex was present in the graph
        """
        removed = False
        complimentary_dictionaries = [(self.parents_map, self.children_map), (self.children_map, self.parents_map)]
        for (first, second) in complimentary_dictionaries:
            if vertex in first:
                keys = first[vertex]
                for key in keys:
                    if key in second:
                        second_values_list = second[key]
                        try:
                            second_values_list.remove(vertex)
                        except ValueError:
                            pass
                        if not second_values_list:
                            second.pop(key)
                first.pop(vertex)
                removed = True
        if removed:
            self.vertices_number -= 1
        return removed

    def add_child(self, parent: int, child: int):
        """!
        @brief This function updated the children map for the given parent-child relationship.
        @param parent The parent id.
        @param child The child id.
        """
        self.children_map[parent].append(child)
        self.vertices_number += 1

    @staticmethod
    def get_graph_from_tree(tree: Tree):
        """!
        @brief This function builds a simple graph from the given tskit tree. Every node in the tree is treated as
        haploid.
        """
        graph = SimpleGraph()
        graph.parents_map = tree.parent_dict
        for (child, parent) in tree.parent_dict.items():
            graph.add_child(parent=parent, child=child)
        return graph

    def add_edge(self, parent: int, child: int):
        """!
        @brief This function updated the children map for the given parent-child relationship.
        @param parent The parent id.
        @param child The child id.
        """
        if parent not in self.children_map:
            self.children_map[parent] = list()
        if child not in self.parents_map:
            self.parents_map[child] = list()
        self.children_map[parent].append(child)
        self.parents_map[child].append(parent)

    @staticmethod
    def get_graph_from_file(filename: str, ploidy: int, max_parent_number: int = 2,
                            missing_parent_notation=None, separation_symbol=' ', skip_first_line: bool = False) \
            -> SimpleGraph:
        """!
        @brief Parses the genealogical graph from the file specified by the filename. The ploidy parameter specifies
        whether the organisms in the file should be treated as haploid or diploid.
        Notice that the every line of the file must contain at least ploidy + 1 ids, but it can optionally have
        some metadata which is ignored by this class.
        @param filename The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param ploidy: The number of ploids that an organism possesses. Must be either 1 or 2.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @return The processed pedigree.

        """
        if ploidy != 1 and ploidy != 2:
            raise Exception(f"The ploidy must be either 1 or 2, found {ploidy} specified")
        pedigree: SimpleGraph = SimpleGraph()

        def process_line(file_line: str):
            if ploidy == 1:
                pedigree.add_haploid_line(line=file_line, max_parent_number=max_parent_number,
                                          missing_parent_notation=missing_parent_notation,
                                          separation_symbol=separation_symbol)
            else:
                pedigree.add_line_from_pedigree(line=file_line,
                                                max_parent_number=max_parent_number,
                                                missing_parent_notation=missing_parent_notation,
                                                separation_symbol=separation_symbol)

        file = open(filename, 'r')
        lines = file.readlines()
        if skip_first_line or lines[0].__contains__('#'):
            lines.pop(0)
        for line in lines:
            process_line(line)
        file.close()
        return pedigree

    @staticmethod
    def get_haploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> SimpleGraph:
        """!
        @brief This method processes the input graph considering that every individual is diploid.
        """
        return SimpleGraph.get_graph_from_file(filename=filename, ploidy=1,
                                               missing_parent_notation=missing_parent_notation,
                                               separation_symbol=separation_symbol,
                                               skip_first_line=skip_first_line)

    @staticmethod
    def get_diploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> SimpleGraph:
        """!
        @brief Parses the pedigree from the file specified by the filename. Every individual is treated as a diploid
        organism.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param filename The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @return The processed pedigree.
        """
        return SimpleGraph.get_graph_from_file(filename=filename, ploidy=2,
                                               max_parent_number=max_parent_number,
                                               missing_parent_notation=missing_parent_notation,
                                               separation_symbol=separation_symbol,
                                               skip_first_line=skip_first_line)

    @staticmethod
    def parse_line(line: str, max_parent_number: int, missing_parent_notation: [str], separation_symbol=' '):
        return list(map(lambda name: int(name) if name not in missing_parent_notation else name,
                        line.strip('\n').split(separation_symbol)[:max_parent_number + 1]))

    def add_haploid_line(self, line: str, max_parent_number: int, separation_symbol=' ', missing_parent_notation=None):
        """!
        @brief Processes the given line and updated the graph, treating every individual as haploid.
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        try:
            child, *parents = SimpleGraph.parse_line(line=line,
                                                     missing_parent_notation=missing_parent_notation,
                                                     max_parent_number=max_parent_number,
                                                     separation_symbol=separation_symbol)
        except ValueError:
            raise Exception("Invalid line")
        if child in self.parents_map:
            warnings.warn(f"Individual {child} is specified multiple times in the graph."
                          f"The previous parents are {self.parents_map[child]}, new values: {parents}", UserWarning)
        for parent in parents:
            if parent not in missing_parent_notation:
                self.add_edge(parent=parent, child=child)

    def add_line_from_pedigree(self, line: str, max_parent_number: int,
                               missing_parent_notation=None, separation_symbol=' '):
        """!
        @brief This function processes a single line from a pedigree and updates the graph accordingly.
        It treats every id as a diploid individual, so if the individual's id is x, then the resulting graph
        will have two vertices 2 * x and 2 * x + 1 for their ploids. If you want to find the individual id by its ploid,
        you can use integer division and divide the id by 2.
        For example, if the given line is "1 2 3" (representing that the individual with id 1 has two parents with
        ids 2 and 3), then the resulting graph will have the following information:
        parents_map[2] = 4, 5
        parents_map[3] = 6, 7
        @param line The line to be parsed. The line must consists of at least three integer values separated by
        the separation symbol.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @param separation_symbol: The symbol used to separate the integers in the line. By default, a space is used.
        @param max_parent_number The maximum number of parents a vertex can have. Must be either 1 or 2
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        child, *parents = SimpleGraph.parse_line(line=line,
                                                 missing_parent_notation=missing_parent_notation,
                                                 max_parent_number=max_parent_number,
                                                 separation_symbol=separation_symbol)
        child_ploid = 2 * int(child)
        for parent in parents:
            if parent not in missing_parent_notation:
                if child_ploid in self.parents_map:
                    raise ValueError("The same individual is specified multiple times in the input file")
                parent = int(parent)
                self.add_edge(2 * parent, child_ploid)
                self.add_edge(2 * parent + 1, child_ploid)
                child_ploid += 1

    def get_vertices(self) -> set[int]:
        return set(self.parents_map).union(self.children_map)

    def get_sink_vertices(self):
        """!
        @brief Returns the sink vertices in the graph (that is, the individuals that have parents,
        but don't have children).
        """
        return [x for x in self.parents_map.keys() if x not in self.children_map.keys()]

    def get_orphans(self):
        """!
        @brief Returns the vertices that don't have parents specified.
        """
        return [x for x in self.children_map if x not in self.parents_map]
