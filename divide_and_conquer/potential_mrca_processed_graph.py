"""!
@file potential_mrca_processed_graph.py
@brief Brief description here
"""

import itertools

import networkx

from graph import Graph
from networkx import flow

print_enabled = False


class PotentialMrcaProcessedGraph(Graph):
    """
    This class represents a preprocessed pedigree graph. Apart from having the usual parents and children mappings,
    it also contains the following information about the pedigree:
        1) The assignments of every vertex to its level. The vertex level is defined to be the length of longest
        path from a proband vertex to the vertex being considered.
        2) A dictionary of dictionaries (a matrix) that maps a vertex u to a dictionary whose keys the u's ancestors
        and the values are the "access vertices" through which u can reach the particular ancestor.
    """

    def __init__(self, pedigree: Graph, probands: [int] = None):
        super().__init__(children_map=pedigree.children_map, parents_map=pedigree.parents_map,
                         vertices_number=pedigree.vertices_number, sink_vertices=pedigree.sink_vertices)
        self.pedigree = pedigree
        self.probands = probands
        self.levels = list()
        self.vertex_to_level_map = dict()
        self.vertex_to_ancestor_map: {int: {int: [int]}} = dict()
        self.initialize_vertex_to_level_map()
        self.initialize_potential_mrca_map()
        # self.inference_cache = dict()

    def initialize_vertex_to_level_map(self):
        if self.probands is None:
            self.probands = {x for x in self.parents_map if x not in self.children_map}
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
            self.vertex_to_ancestor_map[vertex] = dict()

    def initialize_potential_mrca_map(self):
        """!
        @brief This method preprocesses the pedigree and stores the so called "access-to-ancestor" matrix which is
        used during the alignment process. This matrix maps a pair (descendant, ancestor) to the list of the ancestor's
        children vertices through which the descendant can reach the ancestor.
         """
        ancestor_to_vertex_map = {vertex: dict() for vertex in self.vertex_to_level_map.keys()}
        tuple_reuse_map = dict()

        def append_child(child_value: int, children_tuple=None):
            """!
            @brief This is a helper function that appends the child_value to the passed children_tuple.
            The main problem within this preprocessing is that different descendants of the same ancestor vertex
            can climb to this ancestor through the same children. Therefore, in this case, we should reuse the same
            children tuple (through which those descendants climb to the ancestor)
            to save a significant amount of memory.
            @param child_value: A new integer to be appended to the tuple.
            @param children_tuple: The tuple to which the value should be appended.
            @return: The resulting tuple where child_value is added as the last element of the new tuple.
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

    def verify_pmrca_for_vertex_pair(self, first_vertex: int, second_vertex: int, common_ancestor: int):
        """!
        @brief This function verifies that the given vertex is a potential mrca for the given pair of vertices.
        @param first_vertex The first vertex from the pair.
        @param second_vertex The second vertex from the pair.
        @param common_ancestor A common ancestor of the given pair of vertices that is being verified.
        @return True if the given vertex is potential mrca for the given pair of vertices, False otherwise.
        """
        return not (len(self.vertex_to_ancestor_map[first_vertex][common_ancestor]) == 1 and
                    self.vertex_to_ancestor_map[first_vertex][common_ancestor] ==
                    self.vertex_to_ancestor_map[second_vertex][common_ancestor])

    def triplet_condition_holds(self, first_vertex: int, second_vertex: int, third_vertex: int,
                                potential_mrca_candidate: int):
        """!
        @brief This function verifies that the potential mrca candidate satisfies all the conditions involving the
        third vertex. More specifically, it verifies that the
        <a href="https://en.wikipedia.org/wiki/Hall%27s_marriage_theorem">Hall's condition</a> is satisfied for the
        triplet (first_vertex, second_vertex, third_vertex).
        @param first_vertex: The first vertex from the triplet.
        @param second_vertex: The second vertex from the triplet.
        @param third_vertex: The third vertex from the triplet.
        @param potential_mrca_candidate: A common ancestor of the given triplet of vertices that satisfies the
        Hall's condition for the (first_vertex, second_vertex) pair.
        @return True if potential mrca candidate is indeed a potential mrca, False otherwise.
        """
        first_vertex_access = self.vertex_to_ancestor_map[first_vertex][potential_mrca_candidate]
        second_vertex_access = self.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
        current_access = set(first_vertex_access).union(second_vertex_access)
        if len(current_access) > 2:
            return True
        third_vertex_access = self.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
        for access_vertex in third_vertex_access:
            if access_vertex not in current_access:
                return True
        return False

    def extend_pair_to_triple(self, first_vertex: int, second_vertex: int, third_vertex: int,
                              potential_mrca_candidate: int):
        """!
        @brief This function verifies that the Hall condition holds for all the sets involving the third vertex.
        Before verifying the actual constraint, this function tries to check for the most common scenarios and avoid
        unnecessary set calculations.
        @param first_vertex The first vertex from the triplet.
        @param second_vertex The second vertex from the triplet.
        @param third_vertex The third vertex from the triplet.
        @param potential_mrca_candidate A common ancestor of the given triplet of vertices that satisfies the
        Hall's condition for the (first_vertex, second_vertex) pair.
        @return True if potential mrca candidate is indeed a potential mrca, False otherwise.
        """
        first_vertex_access = self.vertex_to_ancestor_map[first_vertex][potential_mrca_candidate]
        second_vertex_access = self.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
        third_vertex_access = self.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
        if first_vertex_access == second_vertex_access == third_vertex_access:
            return False
        if len(third_vertex_access) == 1:
            if third_vertex_access == first_vertex_access:
                return False
            if third_vertex_access == second_vertex_access:
                return False
        return self.triplet_condition_holds(first_vertex, second_vertex, third_vertex, potential_mrca_candidate)

    def verify_additional_constraints(self, partially_assigned_vertices: [int], next_vertex: int, potential_mrca: int):
        """!
        @brief
        @param partially_assigned_vertices:
        @param next_vertex:
        @param potential_mrca:
        @return
        """
        if len(partially_assigned_vertices) == 2:
            [first_vertex, second_vertex] = partially_assigned_vertices
            return self.extend_pair_to_triple(first_vertex, second_vertex, next_vertex, potential_mrca)
        raise Exception("Unimplemented")

    def apply_additional_constraints_for_partial_result(self, coalescent_vertex_to_candidates: {int: [int]},
                                                        partially_assigned_vertices: [int], potential_mrcas: [int],
                                                        next_vertex: int,
                                                        ):
        """!
        @brief
        @param coalescent_vertex_to_candidates:
        @param partially_assigned_vertices:
        @param potential_mrcas:
        @param next_vertex:
        @return:
        """
        next_vertex_candidates = coalescent_vertex_to_candidates[next_vertex]
        next_vertex_result = []
        for candidate in next_vertex_candidates:
            candidate_ancestors = set(self.vertex_to_ancestor_map[candidate])
            candidate_result = []
            for potential_mrca in potential_mrcas:
                if (potential_mrca not in candidate_ancestors or
                        not self.verify_additional_constraints(partially_assigned_vertices, next_vertex,
                                                               potential_mrca)):
                    continue
                candidate_result.append(potential_mrca)
            if candidate_result:
                next_partially_assigned_vertices = list(partially_assigned_vertices)
                next_partially_assigned_vertices.append(candidate)
                next_vertex_result.append((next_partially_assigned_vertices, candidate_result))
        return next_vertex_result

    def get_single_vertex_subsets(self, vertex: int, vertex_ancestors):
        """!
        @brief A helper function returning the partial sub-assignment result for the specified vertex with the
        given ancestors.
        @param vertex
        @param vertex_ancestors
        @return The sub-assignment tuple for a single vertex used by the alignment algorithm
        """
        return [(ancestor, [(1, set(self.vertex_to_ancestor_map[vertex][ancestor]))]) for
                ancestor in vertex_ancestors]

    def get_single_vertex_verified_subset_tuple(self, vertex: int, vertex_ancestors):
        """!
        @brief
        @param vertex:
        @param vertex_ancestors:
        @return:
        """
        return (vertex,), self.get_single_vertex_subsets(vertex, vertex_ancestors)

    def filter_common_ancestors_for_vertex_pair(self, first_candidate, second_candidate):
        """!
        @brief
        @param first_candidate:
        @param second_candidate:
        @return:
        """
        verified_ancestors = []
        if first_candidate not in self.parents_map or second_candidate not in self.parents_map:
            return verified_ancestors
        while self.parents_map[first_candidate] == self.parents_map[second_candidate]:
            verified_ancestors.extend(self.parents_map[first_candidate])
            [first_candidate, second_candidate] = self.parents_map[first_candidate]
            if first_candidate not in self.parents_map or second_candidate not in self.parents_map:
                return verified_ancestors
        first_ancestors = self.vertex_to_ancestor_map[first_candidate]
        second_ancestors = self.vertex_to_ancestor_map[second_candidate]
        if len(second_ancestors) > len(first_ancestors):
            first_ancestors, second_ancestors = second_ancestors, first_ancestors
        for ancestor in first_ancestors:
            if (ancestor in second_ancestors and
                    self.verify_pmrca_for_vertex_pair(first_candidate, second_candidate,
                                                      ancestor)):
                verified_ancestors.append(ancestor)
        return verified_ancestors

    def get_pmracs_for_vertex_pair(self, first: int, second: int, coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first:
        @param second:
        @param coalescent_vertex_to_candidates:
        @return
        """
        result = []
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        if len(first_vertex_candidates) > len(second_vertex_candidates):
            first_vertex_candidates, second_vertex_candidates = second_vertex_candidates, first_vertex_candidates
        for first_vertex_candidate in first_vertex_candidates:
            for second_vertex_candidate in second_vertex_candidates:
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(
                    first_vertex_candidate,
                    second_vertex_candidate)
                if verified_ancestors:
                    result.append(((first_vertex_candidate, second_vertex_candidate), verified_ancestors))
        return result

    def get_pmracs_for_vertex_triple_dfs(self, first: int, second: int, third: int,
                                         coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first:
        @param second:
        @param third:
        @param coalescent_vertex_to_candidates:
        @return:
        """
        result = []
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        third_vertex_candidates = coalescent_vertex_to_candidates[third]
        for first_vertex_candidate in first_vertex_candidates:
            first_vertex_ancestors = self.vertex_to_ancestor_map[first_vertex_candidate]
            if not first_vertex_candidate:
                continue
            for second_vertex_candidate in second_vertex_candidates:
                second_vertex_ancestors = self.vertex_to_ancestor_map[second_vertex_candidate]
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(first_vertex_candidate,
                                                                                  second_vertex_candidate)
                if not verified_ancestors:
                    continue
                for third_vertex_candidate in third_vertex_candidates:
                    third_vertex_ancestors = self.vertex_to_ancestor_map[third_vertex_candidate]
                    if not third_vertex_ancestors:
                        continue
                    if (self.parents_map[first_vertex_candidate] == self.parents_map[second_vertex_candidate] ==
                            self.parents_map[third_vertex_candidate]):
                        fully_verified_ancestors = self.parents_map[first_vertex_candidate]
                    else:
                        fully_verified_ancestors = []
                        for verified_ancestor in verified_ancestors:
                            if (verified_ancestor in third_vertex_ancestors and
                                    self.extend_pair_to_triple(first_vertex_candidate, second_vertex_candidate,
                                                               third_vertex_candidate, verified_ancestor)):
                                fully_verified_ancestors.append(verified_ancestor)
                    if fully_verified_ancestors:
                        result.append(((first_vertex_candidate, second_vertex_candidate, third_vertex_candidate),
                                       fully_verified_ancestors
                                       ))
        return result

    def get_pmracs_for_vertex_triple_iterative(self, first: int, second: int, third: int,
                                               coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief
        @param first:
        @param second:
        @param third:
        @param coalescent_vertex_to_candidates:
        @return:
        """
        first_second_pair_result = self.get_pmracs_for_vertex_pair(first, second, coalescent_vertex_to_candidates)
        first_third_pair_result = self.get_pmracs_for_vertex_pair(first, third, coalescent_vertex_to_candidates)
        first_second_dict = dict()
        first_third_dict = dict()
        for (pair_result, dictionary) in ((first_second_pair_result, first_second_dict),
                                          (first_third_pair_result, first_third_dict)):
            for valid_assignment in pair_result:
                ((first_candidate, _), verified_ancestors) = valid_assignment
                if first_candidate not in dictionary:
                    dictionary[first_candidate] = []
                dictionary[first_candidate].append(valid_assignment)
        result = []
        for (first_candidate, first_second_partial_results) in first_second_dict.items():
            if first_candidate not in first_third_dict:
                continue
            for ((_, second_candidate), verified_ancestors_second) in first_second_partial_results:
                assert first_candidate == _
                verified_ancestors_second_parents_restricted = None
                second_candidate_ancestors = self.vertex_to_ancestor_map[second_candidate]
                for first_third_partial_result in first_third_dict[first_candidate]:
                    ((__, third_candidate), verified_ancestors_third) = first_third_partial_result
                    assert first_candidate == __
                    resulting_ancestors = verified_ancestors_second
                    # If the second and the third candidates happen to have the same parents, we can restrict
                    # the candidates space that should be searched. Specifically, the resulting pmracs
                    # must belong to the intersection of the parents' ancestors sets or be the parents themselves
                    if self.parents_map[second_candidate] == self.parents_map[third_candidate]:
                        if self.parents_map[first_candidate] == self.parents_map[second_candidate]:
                            result.append(((first_candidate, second_candidate, third_candidate),
                                           self.parents_map[first_candidate]))
                            continue
                        if verified_ancestors_second_parents_restricted is None:
                            [first_parent, second_parent] = self.parents_map[second_candidate]
                            parents_set = {first_parent, second_parent}
                            first_parent_ancestors = self.vertex_to_ancestor_map[first_parent]
                            second_parent_ancestors = self.vertex_to_ancestor_map[second_parent]
                            parents_common_ancestors = set(first_parent_ancestors).intersection(
                                set(second_parent_ancestors))
                            verified_ancestors_second_parents_restricted = (
                                (parents_set.union(parents_common_ancestors)).intersection(verified_ancestors_second))
                        resulting_ancestors = verified_ancestors_second_parents_restricted
                    # The resulting pmracs should belong to the both partial pmracs sets
                    resulting_ancestors = [x for x in verified_ancestors_third
                                           if x in resulting_ancestors and
                                           self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, x) and
                                           self.triplet_condition_holds(first_candidate, second_candidate,
                                                                        third_candidate, x)]
                    if resulting_ancestors:
                        result.append(((first_candidate, second_candidate, third_candidate), resulting_ancestors))
        return result

    def get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
                                coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief This function calculates the potential most recent common ancestor (pmrca) for the given vertices.
        TODO: Give the definition of a pmrca taking into account various optimizations that have been implemented
        @param vertices_coalescent_ids The ids for the coalescent vertices.
        @param coalescent_vertex_to_candidates: A dictionary mapping an id of a coalescent vertex to the list of
        ids of pedigree vertices that can be assigned to this coalescent vertex.
        @return All the valid pmracs for the given vertices.
        """
        vertices_length = len(vertices_coalescent_ids)
        vertices_coalescent_ids = sorted(vertices_coalescent_ids,
                                         key=lambda child: len(
                                             coalescent_vertex_to_candidates[child]), reverse=False)
        if vertices_length == 2:
            result = self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                     coalescent_vertex_to_candidates)
        elif vertices_length == 3:
            result = self.get_pmracs_for_vertex_triple_iterative(vertices_coalescent_ids[0],
                                                                 vertices_coalescent_ids[1],
                                                                 vertices_coalescent_ids[2],
                                                                 coalescent_vertex_to_candidates)
        else:
            result = self.get_pmracs_for_vertices_with_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        # for assignment in result:
        #     (assigned_children, verified_candidates_list) = assignment
        #     self.inference_cache[assigned_children] = verified_candidates_list
        if not result:
            raise Exception("Failed")
        return vertices_coalescent_ids, result

    # Returns all the potential common ancestors for the specified vertices
    # Note that not all the potential common ancestors are MRCAs
    def get_pmracs_for_vertices_with_memory(self, vertices_coalescent_ids: [int],
                                            coalescent_vertex_to_candidates: {int: [int]}):
        """!
        @brief Calculates the pmracs for the given vertices in an iterative manner storing the calculations done
        during the verification of the previous sub-assignment.
        """
        vertices_length = len(vertices_coalescent_ids)
        # partial_result elements are lists whose elements
        # have the format: (partial_assignment: tuple, [(candidate, [(preimage_length, image)])])
        alignment_result = []
        first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
        partial_result = []
        for first_candidate in first_candidates:
            # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
            # in the coalescent tree for which we are making the inference
            first_candidate_ancestors = [x for x in self.get_vertex_ancestors(first_candidate)
                                         if len(self.children_map[x]) >= vertices_length]
            if first_candidate_ancestors:
                partial_result.append(self.get_single_vertex_verified_subset_tuple(first_candidate,
                                                                                   first_candidate_ancestors))
        for next_coalescent_id in vertices_coalescent_ids[1:]:
            next_result = []
            if print_enabled:
                print(f"Partial result: {len(partial_result)}")
                print(f"Candidates for the next vertex: {len(coalescent_vertex_to_candidates[next_coalescent_id])}")
                print(
                    f"The whole search space is {len(partial_result) * len(coalescent_vertex_to_candidates[next_coalescent_id])}")
            for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
                next_candidate_ancestors = self.get_vertex_ancestors(next_candidate)
                for verified_children_partial_assignment in partial_result:
                    (assigned_children, verified_candidates_list) = verified_children_partial_assignment
                    next_verified_candidates_list = []
                    extended_assigned_children = assigned_children + (next_candidate,)
                    for verified_candidate_tuple in verified_candidates_list:
                        (verified_candidate, verified_subsets) = verified_candidate_tuple
                        if verified_candidate not in next_candidate_ancestors:
                            continue
                        is_correct_assignment = True
                        verified_candidate_access = self.vertex_to_ancestor_map[next_candidate][verified_candidate]
                        new_subsets = []
                        for verified_subset_tuple in verified_subsets:
                            (preimage_size, image) = verified_subset_tuple
                            image: set
                            new_image = image.union(verified_candidate_access)
                            if len(new_image) < preimage_size + 1:
                                is_correct_assignment = False
                                break
                            new_subsets.append((preimage_size + 1, new_image))
                        if not is_correct_assignment:
                            continue
                        new_subsets.extend(verified_subsets)
                        new_subsets.append(
                            (1, set(self.vertex_to_ancestor_map[next_candidate][verified_candidate])))
                        next_verified_candidates_list.append((verified_candidate, new_subsets))
                        # next_result.append((extended_assigned_children, new_subsets))
                    if next_verified_candidates_list:
                        next_result.append((extended_assigned_children, next_verified_candidates_list))
            if print_enabled:
                print(f"The resulting number of correct assignments is: {len(next_result)}")
            partial_result = next_result
        alignment_result.extend(partial_result)
        result = []
        for assignment in alignment_result:
            (assigned_children, verified_candidates_list) = assignment
            candidates = [verified_candidate_tuple[0] for verified_candidate_tuple in verified_candidates_list]
            result.append((assigned_children, candidates))
        if not result:
            raise Exception("No valid assignments found")
        return result

    def verify_mrca_for_vertices(self, mrca: int, descendant_vertices: [int]):
        # A speed-up for the case with two vertices:
        if len(descendant_vertices) == 2 and \
                self.vertex_to_ancestor_map[descendant_vertices[0]][mrca] != \
                self.vertex_to_ancestor_map[descendant_vertices[1]][mrca]:
            return True
        network_graph = networkx.DiGraph()
        source_vertex_label = "s"
        for descendant in descendant_vertices:
            network_graph.add_edge(source_vertex_label, descendant, capacity=1)
            self.add_edges_to_mrca_from_descendants(network_graph, mrca, descendant_vertices)
        return networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=mrca) == len(
            descendant_vertices)

    def add_edges_to_mrca_from_descendants(self, flow_network: networkx.DiGraph, mrca: int, descendant_vertices: [int]):
        """
        :param flow_network:
        :param mrca:
        :param descendant_vertices:
        :return:
        """
        # TODO: Extract the whole inference logic into the graph matcher class to achieve correct separation of concerns
        mrca_access_label = f"{mrca}'"
        flow_network.add_edge(mrca_access_label, mrca, capacity=len(descendant_vertices))
        for descendant in descendant_vertices:
            # descendant_mrca_access_vertices = self.vertex_to_ancestor_map[descendant][mrca]
            # flow_network.add_edges_from([(x, mrca_access_label) for x in descendant_mrca_access_vertices], capacity=1)
            # current_level_vertices = descendant_mrca_access_vertices
            current_level_vertices = [mrca]
            while current_level_vertices:
                next_level_vertices = set()
                for vertex in current_level_vertices:
                    if vertex == descendant:
                        continue
                    neighbors = self.vertex_to_ancestor_map[descendant][vertex]
                    for neighbor in neighbors:
                        vertex_access_label = f"{vertex}'"
                        # We can potentially override the edge's capacity if vertex is the mrca in the
                        # coalescent tree
                        if not flow_network.has_edge(vertex_access_label, vertex):
                            flow_network.add_edge(vertex_access_label, vertex, capacity=1)
                        flow_network.add_edge(neighbor, vertex_access_label, capacity=1)
                    next_level_vertices.update(neighbors)
                current_level_vertices = next_level_vertices

    def get_vertex_ancestors(self, vertex: int):
        """!
        @brief Returns the vertex's ancestors.
        """
        return self.vertex_to_ancestor_map[vertex].keys()
