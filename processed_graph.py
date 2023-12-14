from abc import abstractmethod

from genealogical_graph import GenealogicalGraph


class ProcessedGraph(GenealogicalGraph):

    def __init__(self, genealogical_graph: GenealogicalGraph):
        super().__init__(genealogical_graph)
        self.descendant_writer = genealogical_graph.descendant_writer

    @abstractmethod
    def initialize_vertex_common_ancestry_entry(self, vertex: int):
        pass

    @abstractmethod
    def update_vertex_common_ancestry_entry(self, vertex: int, ancestor_partner: int, parent: int):
        pass

    @abstractmethod
    def is_correct(self):
        pass

    @abstractmethod
    def get_common_ancestors_map(self, vertex: int):
        pass

    def perform_clean_up(self):
        pass

    @abstractmethod
    def notify_processed_level(self, level_index: int):
        pass

    def initialize_vertex_common_ancestry_entry_given_parents(self, child: int, parent_left: int, parent_right: int):
        self.initialize_vertex_common_ancestry_entry(child)
        self.inherit_ancestry_from_parent(child, parent_left)
        self.inherit_ancestry_from_parent(child, parent_right)

    def inherit_ancestry_from_parent(self, vertex: int, parent: int):
        excluded_descendants = self.get_vertex_excluded_descendants(vertex_id=parent, excluded_vertex=vertex)
        for key, value in self.get_common_ancestors_map(parent).items():
            if key in self.get_vertex_descendants(vertex):
                if key in excluded_descendants:
                    self.update_vertex_common_ancestry_entry(vertex, key, parent)
                else:
                    continue
            self.update_vertex_common_ancestry_entry(vertex, key, parent)

    def get_graph_from_genealogical_graph(self, genealogical_graph: GenealogicalGraph):
        levels = genealogical_graph.levels
        # Initializing the common-ancestor map for the founders
        current_level_vertices = levels[len(levels) - 1]
        for founder_vertex in current_level_vertices:
            self.initialize_vertex_common_ancestry_entry(founder_vertex)
        # Initializing the common-ancestor map for all the other levels
        # TODO Use a custom iterator or a queue instead of level-by-level assignment
        processed = set()
        for index in range(len(levels) - 1).__reversed__():
            print(f"Level {index}")
            for vertex in levels[index]:
                if vertex in processed:
                    continue
                if vertex not in genealogical_graph.parents_map:
                    self.initialize_vertex_common_ancestry_entry(vertex)
                else:
                    [left_parent, right_parent] = genealogical_graph.parents_map[vertex]
                    self.initialize_vertex_common_ancestry_entry_given_parents(vertex,
                                                                               left_parent,
                                                                               right_parent)
                processed.add(vertex)
            self.notify_processed_level(index)
        self.perform_clean_up()
        # self.is_correct()
        return self
