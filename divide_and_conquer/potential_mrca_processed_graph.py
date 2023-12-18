import itertools

from graph import Graph

print_enabled = True


class PotentialMrcaProcessedGraph(Graph):

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
        self.inference_cache = dict()

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
        ancestor_to_vertex_map = {vertex: dict() for vertex in self.vertex_to_level_map.keys()}
        tuple_reuse_map = dict()

        def append_child(child_value: int, children_tuple=None):
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

    def verify_pmrca_for_vertex_pair(self, first_vertex: int, second_vertex: int, potential_mrca: int):
        return not (len(self.vertex_to_ancestor_map[first_vertex][potential_mrca]) == 1 and
                    self.vertex_to_ancestor_map[first_vertex][potential_mrca] ==
                    self.vertex_to_ancestor_map[second_vertex][potential_mrca])

    def extend_pair_to_triple(self, first_vertex: int, second_vertex: int, third_vertex: int,
                              potential_mrca: int):
        first_vertex_access = self.vertex_to_ancestor_map[first_vertex][potential_mrca]
        second_vertex_access = self.vertex_to_ancestor_map[second_vertex][potential_mrca]
        third_vertex_access = self.vertex_to_ancestor_map[third_vertex][potential_mrca]
        if len(third_vertex_access) == 1:
            if third_vertex_access == first_vertex_access:
                return False
            if third_vertex_access == second_vertex_access:
                return False
        current_access = set(first_vertex_access).union(second_vertex_access)
        if len(current_access) > 2:
            return True
        for access_vertex in third_vertex_access:
            if access_vertex not in current_access:
                return True
        return False

    def verify_additional_constraints(self, partially_assigned_vertices: [int], next_vertex: int, potential_mrca: int):
        if len(partially_assigned_vertices) == 2:
            [first_vertex, second_vertex] = partially_assigned_vertices
            return self.extend_pair_to_triple(first_vertex, second_vertex, next_vertex, potential_mrca)
        raise Exception("Unimplemented")

    def apply_additional_constraints_for_partial_result(self, coalescent_vertex_to_candidates: {int: [int]},
                                                        partially_assigned_vertices: [int], potential_mrcas: [int],
                                                        next_vertex: int,
                                                        ):
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
        return [(ancestor, [(1, set(self.vertex_to_ancestor_map[vertex][ancestor]))]) for
                ancestor in vertex_ancestors]

    def get_single_vertex_verified_subset_tuple(self, vertex: int, vertex_ancestors):
        return (vertex,), self.get_single_vertex_subsets(vertex, vertex_ancestors)

    def filter_common_ancestors_for_vertex_pair(self, first_candidate, second_candidate):
        candidate_tuple = tuple(sorted([first_candidate, second_candidate]))
        if candidate_tuple in self.inference_cache:
            return self.inference_cache[candidate_tuple]
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
        self.inference_cache[candidate_tuple] = verified_ancestors
        return verified_ancestors

    def get_pmracs_for_vertex_pair(self, first: int, second: int, coalescent_vertex_to_candidates: {int: [int]}):
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

    def get_pmracs_for_vertex_triple(self, first: int, second: int, third: int,
                                     coalescent_vertex_to_candidates: {int: [int]}):
        result = []
        first_vertex_candidates = coalescent_vertex_to_candidates[first]
        second_vertex_candidates = coalescent_vertex_to_candidates[second]
        third_vertex_candidates = coalescent_vertex_to_candidates[third]
        for first_vertex_candidate in first_vertex_candidates:
            for second_vertex_candidate in second_vertex_candidates:
                verified_ancestors = self.filter_common_ancestors_for_vertex_pair(first_vertex_candidate,
                                                                                  second_vertex_candidate)
                if not verified_ancestors:
                    continue
                for third_vertex_candidate in third_vertex_candidates:

                    triple_candidate_tuple = tuple(sorted([first_vertex_candidate, second_vertex_candidate,
                                                           third_vertex_candidate]))
                    if triple_candidate_tuple in self.inference_cache:
                        fully_verified_ancestors = self.inference_cache[triple_candidate_tuple]
                    else:
                        third_vertex_ancestors = self.vertex_to_ancestor_map[third_vertex_candidate]
                        fully_verified_ancestors = []
                        for verified_ancestor in verified_ancestors:
                            if (verified_ancestor in third_vertex_ancestors and
                                    self.extend_pair_to_triple(first_vertex_candidate, second_vertex_candidate,
                                                               third_vertex_candidate, verified_ancestor)):
                                fully_verified_ancestors.append(verified_ancestor)
                        self.inference_cache[triple_candidate_tuple] = fully_verified_ancestors
                    if fully_verified_ancestors:
                        result.append(((first_vertex_candidate, second_vertex_candidate, third_vertex_candidate),
                                       fully_verified_ancestors
                                       ))
        return result

    def get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
                                coalescent_vertex_to_candidates: {int: [int]}):
        vertices_length = len(vertices_coalescent_ids)
        if vertices_length == 2:
            result = self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                     coalescent_vertex_to_candidates)
        elif vertices_length == 3:
            result = self.get_pmracs_for_vertex_triple(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
                                                       vertices_coalescent_ids[2], coalescent_vertex_to_candidates)
        else:
            result = self.get_pmracs_for_vertices_with_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
        # for assignment in result:
        #     (assigned_children, verified_candidates_list) = assignment
        #     self.inference_cache[assigned_children] = verified_candidates_list
        return result

    # Returns all the potential common ancestors for the specified vertices
    # Note that not all the potential common ancestors are MRCAs
    def get_pmracs_for_vertices_with_memory(self, vertices_coalescent_ids: [int],
                                            coalescent_vertex_to_candidates: {int: [int]}):
        vertices_length = len(vertices_coalescent_ids)
        # partial_result elements are lists whose elements
        # have the format: (partial_assignment: tuple, [(candidate, [(preimage_length, image)])])
        alignment_result = []
        vertices_coalescent_ids = sorted(vertices_coalescent_ids,
                                         key=lambda child: len(
                                             coalescent_vertex_to_candidates[child].keys()), reverse=True)

        first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
        partial_result = []
        for first_candidate in first_candidates:
            # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
            # in the coalescent tree for which we are making the inference
            first_candidate_ancestors = [x for x in self.get_vertex_ancestors(first_candidate)
                                         if len(self.children_map[x]) >= vertices_length]
            # first_candidate_ancestors = self.get_vertex_ancestors(first_candidate)
            if first_candidate_ancestors:
                partial_result.append(self.get_single_vertex_verified_subset_tuple(first_candidate,
                                                                                   first_candidate_ancestors))
        for next_coalescent_id in vertices_coalescent_ids[1:]:
            next_result = []
            if print_enabled:
                print(f"Partial result: {len(partial_result)}")
                print(f"Candidates for the next vertex: {len(coalescent_vertex_to_candidates[next_coalescent_id])}")
                print(f"The whole search space is {len(partial_result) *
                                                   len(coalescent_vertex_to_candidates[next_coalescent_id])}")
            for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
                next_candidate_ancestors = self.get_vertex_ancestors(next_candidate)
                for verified_children_partial_assignment in partial_result:
                    (assigned_children, verified_candidates_list) = verified_children_partial_assignment
                    # assigned_children_list = list(assigned_children)
                    # assigned_children_list.append(next_candidate)
                    # candidates_tuple = tuple(sorted(assigned_children_list))
                    # if candidates_tuple in self.inference_cache:
                    #     print("HIT")
                    next_verified_candidates_list = []
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
                        extended_assigned_children = assigned_children + (next_candidate,)
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

    def get_vertex_ancestors(self, vertex: int):
        return self.vertex_to_ancestor_map[vertex].keys()
