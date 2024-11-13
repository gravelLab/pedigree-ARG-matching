from collections import defaultdict

from graph.genealogical_graph import GenealogicalGraph, SimpleGraph, DescendantMemoryCache
from tskit import TreeSequence, Tree


class CoalescentTree(GenealogicalGraph):
    """!
    This is a helper class that is responsible for working with coalescent trees. Apart from the functionality
    of the GenealogicalGraph, it calculates the connected components (clades) of the graph.
    """

    def __init__(self, graph: SimpleGraph = None):
        if graph is None:
            graph = SimpleGraph()
        super().__init__(pedigree=graph)
        self.descendant_writer = DescendantMemoryCache()
        self.initialize_genealogical_graph_from_probands()

    def get_vertex_parent(self, child: int):
        return self.parents_map[child][0]

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
        return CoalescentTree(graph=pedigree)

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

    def get_identity_solution(self):
        return {x: x for x in self.vertex_to_level_map}

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

    def unmerge_polytomy(self, parent: int, first_child: int, second_child: int, recalculate_levels: bool = True):
        def get_new_vertex_id():
            return max(self.vertex_to_level_map.keys()) + 1

        parent_children = list(self.children_map.get(parent, []))
        if len(parent_children) < 3:
            raise ValueError(f"The specified parent {parent} does not form a polytomy, "
                             f"the only children are {parent_children}")
        selected_children = (first_child, second_child)
        for child in selected_children:
            if child not in parent_children:
                raise ValueError(f"The specified child {child} is not a child of {parent}")
        new_vertex_id = get_new_vertex_id()
        for child in selected_children:
            self.add_edge(parent=new_vertex_id, child=child, recalculate_levels=False)
            self.remove_edge(parent=parent, child=child, recalculate_levels=False)
        self.add_edge(parent=parent, child=new_vertex_id, recalculate_levels=recalculate_levels)

    def write_levels_to_file(self, file, levels):
        for level in levels:
            for vertex in level:
                if vertex in self.parents_map:
                    parents = self.parents_map[vertex]
                    file.write(f"{vertex} {' '.join(str(parent) for parent in parents)}\n")

    @staticmethod
    def get_coalescent_tree_from_file(filepath: str, max_parent_number: int = 2,
                                      missing_parent_notation=None, separation_symbol=' ',
                                      skip_first_line: bool = False):
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filepath,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        result = CoalescentTree(graph=pedigree)
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
