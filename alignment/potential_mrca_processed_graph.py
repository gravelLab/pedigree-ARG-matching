from pathlib import Path

from lineagekit.core.PloidPedigree import PloidPedigree


class PotentialMrcaProcessedGraph(PloidPedigree):
    """
    This class represents a preprocessed pedigree graph. Apart from having the usual parents and children mappings,
    it also contains the following information about the pedigree:

    1) The assignments of every vertex to its level. The vertex level is defined to be the length of longest
    path from a proband vertex to the vertex being considered.

    2) A dictionary of dictionaries (a matrix) that maps a vertex u to a dictionary whose keys the u's ancestors
    and the values are the "access vertices" through which u can reach the particular ancestor.
    """

    def __init__(self, initialize_ancestor_maps: bool = True, *args, **kwargs):
        super().__init__(args, kwargs)
        self.vertex_to_ancestor_map: {int: {int: [int]}} = None
        if initialize_ancestor_maps:
            self.initialize_potential_mrca_map()

    @classmethod
    def from_pedigree(cls, pedigree: PloidPedigree, initialize_ancestor_maps: bool = True):
        obj = cls.__new__(cls)  # Create instance without calling __init__
        obj.__dict__.update(pedigree.__dict__)  # Copy all attributes
        if initialize_ancestor_maps:
            obj.initialize_potential_mrca_map()
        return obj

    def initialize_potential_mrca_map(self):
        """
        This method preprocesses the pedigree and stores the so called "access-to-ancestor" matrix which is
        used during the alignment process. This matrix maps a pair (descendant, ancestor) to the list of the ancestor's
        children vertices through which the descendant can reach the ancestor.
         """
        self.vertex_to_ancestor_map = {key: dict() for key in self}
        ancestor_to_vertex_map = {vertex: dict() for vertex in self}
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

        levels = self.get_levels()
        for vertex in levels[1]:
            vertex_children = self.get_children(vertex)
            for child in vertex_children:
                ancestor_to_vertex_map[vertex][child] = append_child(child)
        for level in levels[2:]:
            for vertex in level:
                vertex_children = self.get_children(vertex)
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
    def get_processed_graph_from_file(filepath: str | Path, missing_parent_notation=None, separation_symbol=' ',
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

        pedigree: PloidPedigree = PloidPedigree.get_ploid_pedigree_from_file(
            filepath=filepath,
            missing_parent_notation=missing_parent_notation,
            separation_symbol=separation_symbol,
            skip_first_line=skip_first_line,
            probands=probands
        )
        return PotentialMrcaProcessedGraph.from_pedigree(pedigree=pedigree, initialize_ancestor_maps=preprocess_graph)

    def get_minimal_path_length(self, descendant: int, ancestor: int) -> int:
        current_level = {ancestor}
        length = 0
        ancestor_map = self.vertex_to_ancestor_map[descendant]
        visited = set(current_level)
        while current_level:
            length += 1
            next_level = {
                v
                for u in current_level
                for v in ancestor_map.get(u, [])
                if v not in visited
            }
            if descendant in next_level:
                return length
            visited.update(next_level)
            current_level = next_level
        raise ValueError(f"Descendant {descendant} not reachable from ancestor {ancestor}")
