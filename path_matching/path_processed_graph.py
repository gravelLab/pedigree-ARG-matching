from bitstring import BitArray

from genealogical_graph import GenealogicalGraph
from graph import Graph
from processed_graph import ProcessedGraph
from itertools import combinations


class PathAwareGenealogicalGraph(GenealogicalGraph):
    class Path:

        # TODO: Append a bit to the path efficiently

        def __init__(self, proband_vertex: int, path: BitArray = None):
            if path is None:
                path = BitArray()
            self.vertex = proband_vertex
            self.path = path
            self.banned = False

        def diverge_from_path(self, parent: int):
            copy_path: BitArray = BitArray(self.path)
            copy_path.append(BitArray(bin=str(parent % 2)))
            # Who created this?
            result = self.__class__(self.vertex, copy_path)
            return result

        def diverge_from_path_reverse(self, parent: int):
            # copy_path: BitArray = BitArray(bin=str(parent % 2))
            # copy_path.append(self.path)
            alternative = BitArray(self.path)
            alternative.prepend(BitArray(1))
            alternative.set(value=(parent % 2), pos=0)
            # Who created Python?
            return self.__class__(self.vertex, alternative)

        def append_parent(self, parent: int):
            self.path.append(BitArray(bin=str(parent % 2)))
            return self

        def ban(self):
            self.banned = True

        def is_banned(self):
            return self.banned

    def __init__(self, pedigree: Graph = None, levels: [[int]] = None):
        super().__init__(pedigree=pedigree)
        self.path_map = dict()

    def additional_process_level_vertex(self, parent_label: int, child: int, level: int):
        # TODO: Redesign the class hierarchy, so that this can be avoided
        if parent_label not in self.path_map:
            self.path_map[parent_label] = list()
            self.path_map[parent_label].append(PathAwareGenealogicalGraph.Path(parent_label))
        paths = self.path_map[child]
        for path in paths:
            self.path_map[parent_label].append(path.diverge_from_path_reverse(parent_label % 2))

    def additional_process_proband_vertex(self, proband_label: int):
        self.path_map.update({proband_label: [self.Path(proband_label)]})

    def initialize_genealogical_graph_from_probands(self, probands: [int] = None):
        super().initialize_genealogical_graph_from_probands()
        # Building the matrix
        result = dict()
        for vertex in self.path_map.keys():
            # TODO: skip for the probands (small optimization)
            dictionary = dict()
            result.update({vertex: dictionary})
            for path in self.path_map[vertex]:
                if path.vertex not in dictionary:
                    dictionary.update({path.vertex: [path]})
                else:
                    dictionary[path.vertex].append(path)
        self.path_map = result


class PathProcessedGraph(ProcessedGraph):

    class PathCommonAncestry:

        def __init__(self, common_ancestor: int, first_vertex: int, second_vertex: int,
                     first_vertex_paths: [PathAwareGenealogicalGraph.Path],
                     second_vertex_paths: [PathAwareGenealogicalGraph.Path]):
            self.common_ancestor = common_ancestor
            self.first_vertex = first_vertex
            self.second_vertex = second_vertex
            self.first_vertex_paths = first_vertex_paths
            self.second_vertex_paths = second_vertex_paths

        def get_paths_for_vertex(self, vertex: int):
            if vertex == self.first_vertex:
                return self.first_vertex_paths
            if vertex == self.second_vertex:
                return self.second_vertex_paths
            raise Exception("Incorrect vertex label")

        def get_other_vertex(self, vertex: int):
            if vertex == self.first_vertex:
                return self.second_vertex
            if vertex == self.second_vertex:
                return self.first_vertex
            raise Exception("Incorrect vertex label")

        def get_paths_for_other_vertex(self, vertex: int):
            return self.get_paths_for_vertex(self.get_other_vertex(vertex))

        def add_paths_for_vertex(self, vertex: int, paths: [PathAwareGenealogicalGraph.Path]):
            if vertex == self.first_vertex:
                self.first_vertex_paths.extend(paths)
            elif vertex == self.second_vertex:
                self.second_vertex_paths.extend(paths)
            else:
                raise Exception("Incorrect vertex label")

        def add_paths_for_other_vertex(self, vertex: int, paths: [PathAwareGenealogicalGraph.Path]):
            self.add_paths_for_vertex(self.get_other_vertex(vertex), paths)

        def update_paths_for_parent(self, parent_common_ancestry, parent: int):
            other_vertex = parent_common_ancestry.get_other_vertex(parent)
            if other_vertex != self.first_vertex and other_vertex != self.second_vertex:
                raise Exception("Incorrect parent label")
            self.add_paths_for_vertex(other_vertex, parent_common_ancestry.get_paths_for_vertex(other_vertex))
            new_paths = [x.diverge_from_path(parent) for x in parent_common_ancestry.get_paths_for_vertex(parent)]
            self.add_paths_for_other_vertex(other_vertex, new_paths)

        def get_inherited_common_ancestry_from_parent(self, parent: int, child: int):
            other_vertex = self.get_other_vertex(parent)
            child_common_ancestry = self.__class__(common_ancestor=self.common_ancestor,
                                                   first_vertex=child,
                                                   second_vertex=other_vertex,
                                                   first_vertex_paths=[],
                                                   second_vertex_paths=[])
            child_common_ancestry.update_paths_for_parent(self, parent)
            return child_common_ancestry

        @staticmethod
        def create_path_common_ancestry_for_founder_vertex(founder_vertex: int, descendant_vertex: int,
                                                           paths: [PathAwareGenealogicalGraph.Path]):
            return PathProcessedGraph.PathCommonAncestry(
                founder_vertex, founder_vertex,
                descendant_vertex,
                [PathAwareGenealogicalGraph.Path(founder_vertex)], paths)

    def __init__(self, genealogical_graph: PathAwareGenealogicalGraph):
        super().__init__(genealogical_graph)
        self.path_map = genealogical_graph.path_map
        self.path_common_ancestry_map = dict()

    def initialize_vertex_common_ancestry_entry(self, vertex: int):
        if vertex in self.path_common_ancestry_map:
            return
        founder_dict = dict()
        self.path_common_ancestry_map.update({vertex: founder_dict})
        for descendant_label in self.path_map[vertex]:
            founder_label = vertex
            paths = self.path_map[vertex][descendant_label]
            path_common_ancestry = PathProcessedGraph.PathCommonAncestry \
                .create_path_common_ancestry_for_founder_vertex(founder_vertex=founder_label,
                                                                descendant_vertex=descendant_label,
                                                                paths=paths)
            last_level_dictionary = dict({vertex: path_common_ancestry})
            founder_dict[descendant_label] = last_level_dictionary
        # self.path_map.pop(vertex)

    def update_vertex_common_ancestry_entry(self, vertex: int, common_ancestor_partner: int, parent: int):
        # The values in this dictionary are dictionary of dictionaries (a matrix) that maps
        # a pair (b, c) to all the path common ancestry objects for the triple (vertex, b, c)
        child_dictionary = self.path_common_ancestry_map[vertex]
        if common_ancestor_partner not in child_dictionary:
            child_dictionary[common_ancestor_partner] = dict()
        partner_map = child_dictionary[common_ancestor_partner]
        for common_ancestor, common_ancestry in self.path_common_ancestry_map[parent][common_ancestor_partner].items():
            if common_ancestor not in partner_map:
                partner_map.update(
                    {common_ancestor: common_ancestry.get_inherited_common_ancestry_from_parent(
                        parent=parent,
                        child=vertex
                    )})
            else:
                partner_map[common_ancestor]: PathProcessedGraph.PathCommonAncestry
                partner_map[common_ancestor].update_paths_for_parent(common_ancestry, parent)

    def is_correct(self):
        pass

    def notify_processed_level(self, level_index: int):
        # keys = (set(self.path_common_ancestry_map.keys()).
                # difference(set(map(lambda vertex: vertex.label, self.levels[level_index]))))
        # for key in keys:
            # self.path_common_ancestry_map.pop(key)
        # TODO The code has been commented out because of the "vertex belonging to level" conflict.
        # TODO The clean-up step can be done when it has been solved
        pass

    def get_common_ancestries_for_vertices(self, vertices: [int]):
        common_ancestors_lists = list()
        for (a, b) in combinations(vertices, 2):
            common_ancestors = self.path_common_ancestry_map[a][b].keys()
            common_ancestors_lists.append(common_ancestors)
        all_vertices_common_ancestors = set.intersection(*common_ancestors_lists)
        return all_vertices_common_ancestors

    def get_common_ancestors_map(self, vertex: int):
        return self.path_common_ancestry_map[vertex]

    def get_common_ancestors_for_vertex_pair(self, first_vertex: int, second_vertex: int):
        if (first_vertex not in self.path_common_ancestry_map or
                second_vertex not in self.path_common_ancestry_map[first_vertex]):
            return dict()
        return self.path_common_ancestry_map[first_vertex][second_vertex]

    def get_vertex_ancestors(self, vertex: int):
        vertex_ancestors = set()
        current_level_ancestors = {vertex}
        previous_length = -1
        new_length = 0
        while previous_length != new_length:
            next_ancestors = set()
            for ancestor in current_level_ancestors:
                if ancestor in self.parents_map:
                    next_ancestors.update(self.parents_map[ancestor])
            vertex_ancestors.update(next_ancestors)
            current_level_ancestors = next_ancestors
            previous_length = new_length
            new_length = len(vertex_ancestors)
        return vertex_ancestors

    def perform_clean_up(self):
        self.path_common_ancestry_map: {int: dict}
        for key in self.path_common_ancestry_map:
            self.path_common_ancestry_map[key].pop(key)
            for common_ancestor_partner in self.path_common_ancestry_map[key]:
                if key in self.path_common_ancestry_map[key][common_ancestor_partner]:
                    self.path_common_ancestry_map[key][common_ancestor_partner].pop(key)
                    self.path_common_ancestry_map[common_ancestor_partner][key].pop(key)
        pass

    def get_common_ancestors(self, vertices: [int]):
        sets = list()
        fixed_vertex = vertices.pop()
        for vertex in vertices:
            for proband_descendant in self.get_proband_descendants(vertex_id=vertex):
                common_ancestors = set(self.get_common_ancestors_for_vertex_pair(
                    first_vertex=proband_descendant,
                    second_vertex=fixed_vertex).keys())
                sets.append(common_ancestors)
        # TODO: Alternatively, find the smallest set in the list
        result = sets.pop()
        for set_entry in sets:
            result = result.intersection(set_entry)
        return result
