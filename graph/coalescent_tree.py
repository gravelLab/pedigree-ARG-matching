import copy
from collections import defaultdict

from tskit import TreeSequence, Tree

from graph.genealogical_graph import GenealogicalGraph, SimpleGraph

from enum import Enum


class Direction(Enum):
    UP = 1,
    DOWN = 2


class CoalescentTree(GenealogicalGraph):
    """!
    This is a helper class that is responsible for working with coalescent trees. Apart from the functionality
    of the GenealogicalGraph, it calculates the connected components (clades) of the graph.
    """

    def __init__(self, graph: SimpleGraph = None, initialize_levels: bool = True):
        if graph is None:
            graph = SimpleGraph()
        super().__init__(pedigree=graph, initialize_levels=initialize_levels)
        # TODO: Refactor the code, so that the descendants can be calculated without the levels
        self.initialize_vertex_to_level_map()
        self.initialize_genealogical_graph_from_probands()

    def clone(self):
        return copy.deepcopy(self)

    def get_vertex_parent(self, child: int):
        """
        This function returns the unique parent vertex of the given vertex.
        Args:
            child: The child vertex id.

        Returns:
            The parent vertex id.
        """
        vertex_parents = self.parents_map.get(child, [])
        if len(vertex_parents) > 1:
            raise Exception(f"The tree is invalid, the vertex {child} has multiple parents")
        if not vertex_parents:
            return None
        return vertex_parents[0]

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

    def get_root_vertex(self):
        top_level_vertices = self.get_top_level_vertices()
        if len(top_level_vertices) > 1:
            raise ValueError("The tree consists of more than one clade")
        return top_level_vertices[0]

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
        pedigree = SimpleGraph(children_map=defaultdict(list, children_map_vertices),
                               parents_map=defaultdict(list, parents_map_vertices))
        return CoalescentTree(graph=pedigree)

    def remove_unary_nodes(self):
        """
        Removes all the unary nodes in the coalescent tree and recalculates the levels.
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

    def get_identity_solution(self):
        return {x: x for x in self.vertex_to_level_map}

    def remove_isolated_vertices(self, recalculate_levels: bool = True):
        isolated_vertices = [x for x in set(self.parents_map).union(self.children_map) if
                             not self.parents_map.get(x, []) and
                             not self.children_map.get(x, [])]
        for isolated_vertex in isolated_vertices:
            self.children_map.pop(isolated_vertex, 0)
            self.parents_map.pop(isolated_vertex, 0)
        if recalculate_levels:
            self.initialize_vertex_to_level_map()

    def merge_edge(self, parent: int, child: int, recalculate_levels: bool = True):
        child_children = list(self.children_map.get(child, []))
        if not child_children:
            raise ValueError(f"The specified {child}-{parent} edge is a proband edge")
        for child_child in child_children:
            self.add_edge(parent=parent, child=child_child, recalculate_levels=False)
            self.remove_edge(parent=child, child=child_child, recalculate_levels=False)
        self.remove_edge(parent=parent, child=child, recalculate_levels=recalculate_levels)

    def unmerge_polytomy(self, child: int, recalculate_levels: bool = True) -> int:
        def get_new_vertex_id():
            return max(self.vertex_to_level_map.keys()) + 1

        child_parent = self.parents_map.get(child, [])
        if not child_parent:
            raise Exception(f"The specified vertex {child} does not have a parent")
        child_parent = child_parent[0]
        child_parent_children = self.children_map[child_parent]
        if len(child_parent_children) < 3:
            raise Exception(f"The specified vertex {child} is not a part of a polytomy")
        self.remove_edge(parent=child_parent, child=child, recalculate_levels=False)
        new_vertex_id = get_new_vertex_id()
        self.add_edge(parent=new_vertex_id, child=child, recalculate_levels=False)
        child_parent_parent = self.parents_map.get(child_parent, [])
        if child_parent_parent:
            child_parent_parent = child_parent_parent[0]
            self.parents_map[new_vertex_id] = []
            self.remove_edge(parent=child_parent_parent, child=child_parent, recalculate_levels=False)
            self.add_edge(parent=child_parent_parent, child=new_vertex_id, recalculate_levels=False)
        else:
            self.parents_map[child_parent] = []
        self.add_edge(parent=new_vertex_id, child=child_parent, recalculate_levels=recalculate_levels)
        return new_vertex_id

    def write_levels_to_file(self, file, levels):
        for level in levels:
            for vertex in level:
                vertex_parents = self.parents_map.get(vertex, ())
                if vertex_parents:
                    file.write(f"{vertex} {' '.join(str(parent) for parent in vertex_parents)}\n")

    def cut_edge(self, edge_child_vertex: int, direction: Direction):
        bottom_tree = self.get_vertex_descendants(edge_child_vertex)
        # The chosen direction specifies which part of the tree should be kept
        if direction == Direction.UP:
            vertices_to_remove = bottom_tree
        elif direction == Direction.DOWN:
            vertices_to_remove = self.get_vertices().difference(bottom_tree)
        else:
            raise ValueError(f"The specified direction {direction} is not supported")
        for vertex in vertices_to_remove:
            self.remove_vertex(vertex, recalculate_levels=False)
        self.initialize_vertex_to_level_map()
        self.remove_unary_nodes()

    @staticmethod
    def get_coalescent_tree_from_file(filepath: str,
                                      missing_parent_notation=None, separation_symbol=' ',
                                      skip_first_line: bool = False, initialize_levels: bool = True):
        graph: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filepath,
                                                                     max_parent_number=1,
                                                                     missing_parent_notation=missing_parent_notation,
                                                                     separation_symbol=separation_symbol,
                                                                     skip_first_line=skip_first_line)
        result = CoalescentTree(graph=graph, initialize_levels=initialize_levels)
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
