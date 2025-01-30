from graph.genealogical_graph import GenealogicalGraph, SimpleGraph


class PotentialMrcaProcessedGraph(GenealogicalGraph):
    """
    This class represents a preprocessed pedigree graph. Apart from having the usual parents and children mappings,
    it also contains the following information about the pedigree:

    1) The assignments of every vertex to its level. The vertex level is defined to be the length of longest
    path from a proband vertex to the vertex being considered.

    2) A dictionary of dictionaries (a matrix) that maps a vertex u to a dictionary whose keys the u's ancestors
    and the values are the "access vertices" through which u can reach the particular ancestor.
    """

    def __init__(self, pedigree: SimpleGraph, probands: [int] = None, initialize_ancestor_maps: bool = True,
                 initialize_levels: bool = True):
        super().__init__(pedigree=pedigree, probands=probands, initialize_levels=initialize_levels)
        self.vertex_to_ancestor_map: {int: {int: [int]}} = {key: dict() for key in self.vertex_to_level_map}
        if initialize_ancestor_maps:
            self.initialize_potential_mrca_map()
        # self.inference_cache = dict()

    def initialize_potential_mrca_map(self):
        """
        This method preprocesses the pedigree and stores the so called "access-to-ancestor" matrix which is
        used during the alignment process. This matrix maps a pair (descendant, ancestor) to the list of the ancestor's
        children vertices through which the descendant can reach the ancestor.
         """
        # TODO: Recalculate the levels if necessary
        self.vertex_to_ancestor_map = {key: dict() for key in self.vertex_to_level_map}
        ancestor_to_vertex_map = {vertex: dict() for vertex in self.vertex_to_level_map.keys()}
        tuple_reuse_map = dict()

        def append_child(child_value: int, children_tuple=None):
            """
            This is a helper function that appends the child_value to the passed children_tuple.
            The main problem within this preprocessing is that different descendants of the same ancestor vertex
            can climb to this ancestor through the same children. Therefore, in this case, we should reuse the same
            children tuple (through which those descendants climb to the ancestor)
            to save a significant amount of memory.

            Args:
                child_value: A new integer to be appended to the tuple.
                children_tuple: The tuple to which the value should be appended.

            Returns:
                The resulting tuple where child_value is added as the last element of the new tuple.
            """
            if children_tuple is None:
                children_tuple = tuple()
            new_tuple = children_tuple + (child_value,)
            if new_tuple in tuple_reuse_map:
                return tuple_reuse_map[new_tuple]
            tuple_reuse_map[new_tuple] = new_tuple
            return new_tuple

        for vertex in self.levels[1]:
            vertex_children = self.children_map[vertex]
            for child in vertex_children:
                ancestor_to_vertex_map[vertex][child] = append_child(child)
        for level in self.levels[2:]:
            for vertex in level:
                vertex_children = self.children_map[vertex]
                for child in vertex_children:
                    ancestor_to_vertex_map[vertex][child] = append_child(child)
                    for child_descendant in ancestor_to_vertex_map[child].keys():
                        if child_descendant in ancestor_to_vertex_map[vertex]:
                            ancestor_to_vertex_map[vertex][child_descendant] = (
                                append_child(child, ancestor_to_vertex_map[vertex][child_descendant]))
                        else:
                            ancestor_to_vertex_map[vertex][child_descendant] = append_child(child)
        for ancestor, descendant_map in ancestor_to_vertex_map.items():
            descendant_map: dict
            for descendant, access_vertices in descendant_map.items():
                self.vertex_to_ancestor_map[descendant][ancestor] = access_vertices
            # Perform incremental clean-up
            descendant_map.clear()
        pass

    def get_vertex_ancestors(self, vertex: int):
        """
        Returns the vertex's ancestors.

        Args:
            vertex (int): The vertex for which the ancestors should be returned.
        """
        return self.vertex_to_ancestor_map[vertex].keys()

    @staticmethod
    def get_processed_graph_from_file(filepath: str, missing_parent_notation=None, separation_symbol=' ',
                                      preprocess_graph: bool = True,
                                      probands: [int] = None,
                                      skip_first_line: bool = False
                                      ):
        """
        Parses the pedigree file and builds the graph.

        Args:

            filepath: The path to the pedigree file.
            missing_parent_notation: String representing an unspecified parent.
            separation_symbol: The sting used to separate columns in the input file.
            preprocess_graph: Specifies whether the graph should be preprocessed for the alignment algorithm.
            probands: If specified, the parsed graph is reduced to the ascending genealogy of the specified probands.
            skip_first_line: Specifies whether the first line should be skipped.

        Returns:
        """
        pedigree: SimpleGraph = SimpleGraph.get_diploid_graph_from_file(filepath=filepath,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        if probands is None:
            return PotentialMrcaProcessedGraph(pedigree=pedigree, initialize_levels=preprocess_graph,
                                               initialize_ancestor_maps=preprocess_graph)
        unprocessed_graph = PotentialMrcaProcessedGraph(pedigree=pedigree, initialize_levels=False,
                                                        initialize_ancestor_maps=False)
        unprocessed_graph.reduce_to_ascending_genealogy(probands=probands, recalculate_levels=preprocess_graph)
        if preprocess_graph:
            unprocessed_graph.initialize_potential_mrca_map()
        return unprocessed_graph
