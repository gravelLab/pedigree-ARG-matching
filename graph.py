"""!
@file graph.py
@brief The file contains the realization of the basic Graph class that is responsible for parsing a pedigree file
and building the parents and children maps.
"""


class Graph:
    """!
    This class represents a simple genealogical graph and stores the information about the vertices and the edges
    (children-parent relationships) in a graph
    """

    def __init__(self, children_map=None, parents_map=None, sink_vertices=None, vertices_number=0):
        if parents_map is None:
            parents_map = dict()
        if children_map is None:
            children_map = dict()
        self.children_map = children_map
        self.parents_map = parents_map
        self.vertices_number = vertices_number
        self.sink_vertices = sink_vertices
        # TODO: Add the "vertices" set

    def get_connected_component_for_vertex(self, vertex: int):
        """!
        @brief The function finds all the vertices that are in the same connected component as the passed vertex.
        @param vertex The vertex for which the connected component should be calculated.
        @return The set of vertices that are in the same connected component as the passed vertex.
        """
        result = {vertex}
        previous_length = 0
        new_length = 1
        while previous_length != new_length:
            new_result = set(result)
            for connected_component_vertex in result:
                if connected_component_vertex in self.children_map:
                    new_result.update(self.children_map[connected_component_vertex])
                if connected_component_vertex in self.parents_map:
                    new_result.update(self.parents_map[connected_component_vertex])
            previous_length = new_length
            new_length = len(new_result)
            result = new_result
        return result

    def remove_vertex(self, vertex: int):
        """!
        @brief This function removes the vertex from the graph.
        @param vertex The vertex to be removed.
        """
        removed = False
        if vertex in self.parents_map:
            self.parents_map.pop(vertex)
            removed = True
        if vertex in self.children_map:
            self.children_map.pop(vertex)
            removed = True
        if removed:
            self.vertices_number -= 1
        if self.sink_vertices and vertex in self.sink_vertices:
            self.sink_vertices.remove(vertex)

    def add_child(self, parent: int, child: int):
        """!
        @brief This function updated the children map for the given parent-child relationship.
        @param parent The parent id.
        @param child The child id.
        """
        if parent in self.children_map:
            self.children_map[parent].append(child)
        else:
            self.children_map[parent] = list([child])
        self.vertices_number += 1

    def add_line_from_pedigree(self, line, missing_parent_notation=None, separation_symbol=' '):
        """!
        @brief This function processes a single line from a pedigree and updates the graph accordingly.
        @param line The line to be parsed. The line must consists of at least three integer values separated by
        the separation symbol.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @param separation_symbol: The symbol used to separate the integers in the line. By default, a space is used.
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        (child, mother, father) = list(map(lambda name: str(name), line.strip('\n').split(separation_symbol)))[:3]
        if mother in missing_parent_notation or father in missing_parent_notation:
            return
        (child, mother, father) = (int(child), int(mother), int(father))
        child_maternal = 2 * child
        child_paternal = child_maternal + 1
        self.add_child(2 * mother, child_maternal)
        self.add_child(2 * mother + 1, child_maternal)
        self.add_child(2 * father, child_paternal)
        self.add_child(2 * father + 1, child_paternal)
        # Every non-founder appears only once in the pedigree as a child
        assert child_maternal not in self.parents_map
        assert child_paternal not in self.parents_map
        self.parents_map[child_maternal] = [2 * mother, 2 * mother + 1]
        self.parents_map[child_paternal] = [2 * father, 2 * father + 1]

    def get_sink_vertices(self):
        """!
        @brief Returns the sink vertices in the graph (that is, the individuals that have parents,
        but don't have children).
        """
        if not self.sink_vertices:
            self.sink_vertices = [x for x in self.parents_map.keys() if x not in self.children_map.keys()]
        return self.sink_vertices

    def get_orphans(self):
        """!
        @brief Returns the vertices that don't have parents specified.
        """
        return [x for x in self.children_map.keys() if x not in self.parents_map.keys()]

    @staticmethod
    def get_pedigree_from_file(filename, missing_parent_notation=None, separation_symbol=' '):
        """!
        @brief Parses the pedigree from the file specified by the filename.
        @param filename The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @return The processed pedigree.
        """
        pedigree = Graph()
        file = open(filename, 'r')
        lines = file.readlines()
        if not lines[0].__contains__('#'):
            pedigree.add_line_from_pedigree(lines[0], missing_parent_notation, separation_symbol)
        lines.pop(0)
        for line in lines:
            pedigree.add_line_from_pedigree(line, missing_parent_notation, separation_symbol)
        file.close()
        return pedigree

    # @staticmethod
    # def get_graph_from_tree(tree: Tree):
    #     graph = Graph()
    #     graph.parents_map = tree.parent_dict
    #     for (child, parent) in tree.parent_dict.items():
    #         graph.add_child(parent=parent, child=child)
    #     return graph
