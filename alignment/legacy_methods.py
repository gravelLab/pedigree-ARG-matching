# def get_pmracs_for_vertices_with_maximum_flow_iterative(self, vertices_coalescent_ids: [int],
#                                                         coalescent_vertex_to_candidates: {int: [int]}):
#     vertices_length = len(vertices_coalescent_ids)
#     first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
#     partial_result = []
#     for first_candidate in first_candidates:
#         # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
#         # in the coalescent tree for which we are making the inference
#         first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(first_candidate)
#                                      if len(self.pedigree.children_map[x]) >= vertices_length]
#         if first_candidate_ancestors:
#             partial_result.append((first_candidate, first_candidate_ancestors))
#     for next_coalescent_id in vertices_coalescent_ids[1:]:
#         next_result = []
#         print_log(f"Partial result: {len(partial_result)}")
#         print_log(f"Candidates for the next vertex: {len(coalescent_vertex_to_candidates[next_coalescent_id])}")
#         print_log(f"The whole search space is "
#                   f"{len(partial_result) * len(coalescent_vertex_to_candidates[next_coalescent_id])}")
#         for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
#             next_candidate_ancestors = self.pedigree.get_vertex_ancestors(next_candidate)
#             for verified_children_partial_assignment in partial_result:
#                 (assigned_children, verified_candidates_list) = verified_children_partial_assignment
#                 next_verified_candidates_list = []
#                 extended_assigned_children = assigned_children + (next_candidate,)
#                 for verified_candidate in verified_candidates_list:
#                     if verified_candidate not in next_candidate_ancestors:
#                         continue
#                     if self.verify_mrca_for_vertices(verified_candidate, extended_assigned_children):
#                         next_verified_candidates_list.append(verified_candidate)
#                 if next_verified_candidates_list:
#                     next_result.append((extended_assigned_children, next_verified_candidates_list))
#         print_log(f"The resulting number of correct assignments is: {len(next_result)}")
#         partial_result = next_result
#     if not partial_result:
#         raise Exception("No valid assignments found")
#     return partial_result
#

#     def get_pmrcas_without_memory(self, vertices_coalescent_ids: [int], coalescent_vertex_to_candidates: {int: [int]}):
#         vertices_length = len(vertices_coalescent_ids)
#         first_candidates = coalescent_vertex_to_candidates[vertices_coalescent_ids[0]]
#         partial_result = []
#         for first_candidate in first_candidates:
#             # Filtering out the pedigree candidates who have fewer children than the number of children of the vertex
#             # in the coalescent tree for which we are making the inference
#             first_candidate_ancestors = [x for x in self.pedigree.get_vertex_ancestors(first_candidate)
#                                          if len(self.pedigree.children_map[x]) >= vertices_length]
#             if first_candidate_ancestors:
#                 partial_result.append(([first_candidate], first_candidate_ancestors))
#         for next_coalescent_id in vertices_coalescent_ids[1:]:
#             next_result = []
#             for next_candidate in coalescent_vertex_to_candidates[next_coalescent_id]:
#                 next_candidate_ancestors = self.pedigree.get_vertex_ancestors(next_candidate)
#                 for verified_children_partial_assignment in partial_result:
#                     (assigned_children, verified_candidates_list) = verified_children_partial_assignment
#                     new_assigned_children_list = assigned_children + [next_candidate]
#                     new_verified_candidates_list = []
#                     for candidate in verified_candidates_list:
#                         if candidate not in next_candidate_ancestors:
#                             continue
#                         if self.verify_pmrca_for_vertices(candidate, new_assigned_children_list):
#                             new_verified_candidates_list.append(candidate)
#                     if new_verified_candidates_list:
#                         next_result.append((new_assigned_children_list, new_verified_candidates_list))
#             partial_result = next_result
#         return partial_result
#
#

# def verify_mrca_for_vertices(self, mrca: int, descendant_vertices: [int]):
#     network_graph = networkx.DiGraph()
#     source_vertex_label = "s"
#     for descendant in descendant_vertices:
#         network_graph.add_edge(source_vertex_label, descendant, capacity=1)
#         self.add_edges_to_mrca_from_descendants(network_graph, mrca, descendant_vertices)
#     return networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=mrca) == len(
#         descendant_vertices)
#
#
# def verify_pmrca_for_vertices(self, pmrca: int, descendant_vertices: [int]):
#     # A speed-up for the case with two vertices
#     if len(descendant_vertices) == 2 and \
#             self.pedigree.vertex_to_ancestor_map[descendant_vertices[0]][pmrca] != \
#             self.pedigree.vertex_to_ancestor_map[descendant_vertices[1]][pmrca]:
#         return True
#     network_graph = networkx.DiGraph()
#     source_vertex_label = "s"
#     for descendant in descendant_vertices:
#         network_graph.add_edge(source_vertex_label, descendant, capacity=1)
#         for access_vertex in self.pedigree.vertex_to_ancestor_map[descendant][pmrca]:
#             network_graph.add_edge(descendant, access_vertex, capacity=1)
#             network_graph.add_edge(access_vertex, pmrca, capacity=1)
#     return networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=pmrca) == len(
#         descendant_vertices)

# def get_pmracs_for_vertices(self, vertices_coalescent_ids: [int],
#                             coalescent_vertex_to_candidates: {int: [int]}):
#     """!
#     @brief This function calculates the potential most recent common ancestor (pmrca) for the given vertices.
#     TODO: Give the definition of a pmrca taking into account various optimizations that have been implemented
#     @param vertices_coalescent_ids The ids for the coalescent vertices.
#     @param coalescent_vertex_to_candidates: A dictionary mapping an id of a coalescent vertex to the list of
#     ids of pedigree vertices that can be assigned to this coalescent vertex.
#     @return All the valid pmracs for the given vertices.
#     """
#     vertices_length = len(vertices_coalescent_ids)
#     vertices_coalescent_ids = sorted(vertices_coalescent_ids,
#                                      key=lambda child: len(coalescent_vertex_to_candidates[child]),
#                                      reverse=False
#                                      )
#     if vertices_length == 2:
#         result = self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[0], vertices_coalescent_ids[1],
#                                                  coalescent_vertex_to_candidates)
#     elif vertices_length == 3:
#         result = self.get_pmracs_for_vertex_triple_iterative(vertices_coalescent_ids[0],
#                                                              vertices_coalescent_ids[1],
#                                                              vertices_coalescent_ids[2],
#                                                              coalescent_vertex_to_candidates)
#     else:
#         result = self.get_pmracs_for_vertices_dfs(vertices_coalescent_ids, coalescent_vertex_to_candidates)
#     # else:
#     #     separate_results = []
#     #     start_index = 0
#     #     if len(vertices_coalescent_ids) % 2 == 1:
#     #         separate_results.append(self.get_pmracs_for_vertex_triple_iterative(
#     #             vertices_coalescent_ids[0],
#     #             vertices_coalescent_ids[1],
#     #             vertices_coalescent_ids[2],
#     #             coalescent_vertex_to_candidates
#     #         ))
#     #         start_index = 3
#     #     separate_results.extend([self.get_pmracs_for_vertex_pair(vertices_coalescent_ids[i],
#     #                                                              vertices_coalescent_ids[i + 1],
#     #                                                              coalescent_vertex_to_candidates)
#     #                              for i in range(start_index, len(vertices_coalescent_ids), 2)])
#     #     common_keys = frozenset.intersection(*(frozenset(d.keys()) for d in separate_results))
#     #     separate_results = [{key: d[key] for key in common_keys} for d in separate_results]
#     #     result = separate_results[0]
#     #     for next_result in separate_results[1:]:
#     #         new_result = dict()
#     #         for ancestor_candidate, children_assignments in result.items():
#     #             extended_children_assignments = []
#     #             for children_assignment in children_assignments:
#     #                 for next_children_assignment in next_result[ancestor_candidate]:
#     #                     joined_children_assignment = children_assignment + next_children_assignment
#     #                     if (len(frozenset(next_children_assignment).union(children_assignment))
#     #                             != len(children_assignment) + len(next_children_assignment)):
#     #                         continue
#     #                     if (len(joined_children_assignment) > 3 or
#     #                             self.verify_pmrca_for_vertices(ancestor_candidate, joined_children_assignment)):
#     #                         extended_children_assignments.append(joined_children_assignment)
#     #             if extended_children_assignments:
#     #                 new_result[ancestor_candidate] = extended_children_assignments
#     #         result = new_result
#     # result = self.get_pmrcas_without_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
#     # result = self.get_pmracs_for_vertices_with_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
#     # result = self.get_mrcas_without_memory(vertices_coalescent_ids, coalescent_vertex_to_candidates)
#     # for assignment in result:
#     #     (assigned_children, verified_candidates_list) = assignment
#     #     self.inference_cache[assigned_children] = verified_candidates_list
#     # if not result:
#     #     raise Exception("Failed the alignment, 0 candidates found")
#     # candidates_ancestors = {y for x in result for y in self.pedigree.vertex_to_ancestor_map[x]}
#     # result = {k: v for k, v in result.items() if k not in candidates_ancestors}
#     return vertices_coalescent_ids, result

# def get_pmracs_for_vertex_triple_iterative(self, first: int, second: int, third: int,
#                                            coalescent_vertex_to_candidates: {int: [int]}):
#     """!
#     @brief
#     @param first: The first coalescent vertex
#     @param second:The second coalescent vertex
#     @param third: The third coalescent vertex
#     @param coalescent_vertex_to_candidates: The map returning the candidate list for a given coalescent vertex
#     @return:
#     """
#     triplet_cache = dict()
#
#     def get_triplet_tuple(first_element, second_element, third_element):
#         triplet_tuple = (first_element, second_element, third_element)
#         if triplet_tuple in triplet_cache:
#             return triplet_cache[triplet_tuple]
#         triplet_cache[triplet_tuple] = triplet_tuple
#         return triplet_tuple
#
#     first_second_pair_result = self.get_pmracs_for_vertex_pair(first, second, coalescent_vertex_to_candidates)
#     first_third_pair_result = self.get_pmracs_for_vertex_pair(first, third, coalescent_vertex_to_candidates)
#     first_second_dict = defaultdict(dict)
#     first_third_dict = defaultdict(dict)
#     for (pair_result, dictionary) in ((first_second_pair_result, first_second_dict),
#                                       (first_third_pair_result, first_third_dict)):
#         for valid_assignment in pair_result.items():
#             (pmrca_candidate, children_assignments) = valid_assignment
#             for children_assignment in children_assignments:
#                 (first_candidate, other_candidate) = children_assignment
#                 if first_candidate not in dictionary[pmrca_candidate]:
#                     dictionary[pmrca_candidate][first_candidate] = []
#                 dictionary[pmrca_candidate][first_candidate].append(other_candidate)
#     shared_common_ancestors = first_third_pair_result.keys() & first_third_pair_result.keys()
#     result = defaultdict(list)
#     for shared_common_ancestor in shared_common_ancestors:
#         shared_first_vertex_assignments = (first_second_dict[shared_common_ancestor].keys() &
#                                            first_third_dict[shared_common_ancestor].keys())
#         for first_candidate in shared_first_vertex_assignments:
#             for second_candidate in first_second_dict[shared_common_ancestor][first_candidate]:
#                 for third_candidate in first_third_dict[shared_common_ancestor][first_candidate]:
#                     triplet_tuple = get_triplet_tuple(first_candidate, second_candidate, third_candidate)
#                     # TODO: Add common parents speed-up
#                     if (self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, shared_common_ancestor)
#                             and self.triplet_condition_holds(first_candidate, second_candidate,
#                                                              third_candidate, shared_common_ancestor)):
#                         result[shared_common_ancestor].append(triplet_tuple)
#
#     # result = []
#     # for (first_candidate, first_second_partial_results) in first_second_dict.items():
#     #     if first_candidate not in first_third_dict:
#     #         continue
#     #     for ((_, second_candidate), verified_ancestors_second) in first_second_partial_results:
#     #         assert first_candidate == _
#     #         verified_ancestors_second_parents_restricted = None
#     #         second_candidate_ancestors = self.pedigree.vertex_to_ancestor_map[second_candidate]
#     #         for first_third_partial_result in first_third_dict[first_candidate]:
#     #             ((__, third_candidate), verified_ancestors_third) = first_third_partial_result
#     #             assert first_candidate == __
#     #             resulting_ancestors = verified_ancestors_second
#     #             # If the second and the third candidates happen to have the same parents, we can restrict
#     #             # the candidates space that should be searched. Specifically, the resulting pmracs
#     #             # must belong to the intersection of the parents' ancestors sets or be the parents themselves
#     #             if self.pedigree.parents_map[second_candidate] == self.pedigree.parents_map[third_candidate]:
#     #                 if self.pedigree.parents_map[first_candidate] == self.pedigree.parents_map[second_candidate]:
#     #                     result.append(((first_candidate, second_candidate, third_candidate),
#     #                                    self.pedigree.parents_map[first_candidate]))
#     #                     continue
#     #                 if verified_ancestors_second_parents_restricted is None:
#     #                     [first_parent, second_parent] = self.pedigree.parents_map[second_candidate]
#     #                     parents_set = {first_parent, second_parent}
#     #                     first_parent_ancestors = self.pedigree.vertex_to_ancestor_map[first_parent]
#     #                     second_parent_ancestors = self.pedigree.vertex_to_ancestor_map[second_parent]
#     #                     parents_common_ancestors = set(first_parent_ancestors).intersection(
#     #                         set(second_parent_ancestors))
#     #                     verified_ancestors_second_parents_restricted = (
#     #                         (parents_set.union(parents_common_ancestors)).intersection(verified_ancestors_second))
#     #                 resulting_ancestors = verified_ancestors_second_parents_restricted
#     #             # The resulting pmracs should belong to the both partial pmracs sets
#     #             resulting_ancestors = [x for x in verified_ancestors_third
#     #                                    if x in resulting_ancestors and
#     #                                    self.verify_pmrca_for_vertex_pair(second_candidate, third_candidate, x) and
#     #                                    self.triplet_condition_holds(first_candidate, second_candidate,
#     #                                                                 third_candidate, x)]
#     #             if resulting_ancestors:
#     #                 result.append(((first_candidate, second_candidate, third_candidate), resulting_ancestors))
#     return result

# def get_pmracs_for_vertex_triple_dfs(self, first: int, second: int, third: int,
#                                      coalescent_vertex_to_candidates: {int: [int]}):
#     """!
#     @brief
#     @param first:
#     @param second:
#     @param third:
#     @param coalescent_vertex_to_candidates:
#     @return:
#     """
#     result = []
#     first_vertex_candidates = coalescent_vertex_to_candidates[first]
#     second_vertex_candidates = coalescent_vertex_to_candidates[second]
#     third_vertex_candidates = coalescent_vertex_to_candidates[third]
#     for first_vertex_candidate in first_vertex_candidates:
#         first_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[first_vertex_candidate]
#         if not first_vertex_candidate:
#             continue
#         for second_vertex_candidate in second_vertex_candidates:
#             second_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[second_vertex_candidate]
#             verified_ancestors = self.filter_common_ancestors_for_vertex_pair(first_vertex_candidate,
#                                                                               second_vertex_candidate)
#             if not verified_ancestors:
#                 continue
#             for third_vertex_candidate in third_vertex_candidates:
#                 third_vertex_ancestors = self.pedigree.vertex_to_ancestor_map[third_vertex_candidate]
#                 if not third_vertex_ancestors:
#                     continue
#                 if (self.pedigree.parents_map[first_vertex_candidate] ==
#                         self.pedigree.parents_map[second_vertex_candidate] ==
#                         self.pedigree.parents_map[third_vertex_candidate]):
#                     fully_verified_ancestors = self.pedigree.parents_map[first_vertex_candidate]
#                 else:
#                     fully_verified_ancestors = []
#                     for verified_ancestor in verified_ancestors:
#                         if (verified_ancestor in third_vertex_ancestors and
#                                 self.extend_pair_to_triple(first_vertex_candidate, second_vertex_candidate,
#                                                            third_vertex_candidate, verified_ancestor)):
#                             fully_verified_ancestors.append(verified_ancestor)
#                 if fully_verified_ancestors:
#                     result.append(((first_vertex_candidate, second_vertex_candidate, third_vertex_candidate),
#                                    fully_verified_ancestors
#                                    ))
#     return result


# def get_single_vertex_verified_subset_tuple(self, vertex: int, vertex_ancestors):
#     """!
#     @brief
#     @param vertex:
#     @param vertex_ancestors:
#     @return:
#     """
#     return (vertex,), self.get_single_vertex_subsets(vertex, vertex_ancestors)

# def get_single_vertex_subsets(self, vertex: int, vertex_ancestors):
#     """!
#     @brief A helper function returning the partial sub-assignment result for the specified vertex with the
#     given ancestors.
#     @param vertex
#     @param vertex_ancestors
#     @return The sub-assignment tuple for a single vertex used by the alignment algorithm
#     """
#     return [(ancestor, [(1, frozenset(self.pedigree.vertex_to_ancestor_map[vertex][ancestor]))]) for
#             ancestor in vertex_ancestors]

# def apply_additional_constraints_for_partial_result(self, coalescent_vertex_to_candidates: {int: [int]},
#                                                     partially_assigned_vertices: [int], potential_mrcas: [int],
#                                                     next_vertex: int,
#                                                     ):
#     """!
#     @brief
#     @param coalescent_vertex_to_candidates:
#     @param partially_assigned_vertices:
#     @param potential_mrcas:
#     @param next_vertex:
#     @return:
#     """
#     next_vertex_candidates = coalescent_vertex_to_candidates[next_vertex]
#     next_vertex_result = []
#     for candidate in next_vertex_candidates:
#         candidate_ancestors = set(self.pedigree.vertex_to_ancestor_map[candidate])
#         candidate_result = []
#         for potential_mrca in potential_mrcas:
#             if (potential_mrca not in candidate_ancestors or
#                     not self.verify_additional_constraints(partially_assigned_vertices, next_vertex,
#                                                            potential_mrca)):
#                 continue
#             candidate_result.append(potential_mrca)
#         if candidate_result:
#             next_partially_assigned_vertices = list(partially_assigned_vertices)
#             next_partially_assigned_vertices.append(candidate)
#             next_vertex_result.append((next_partially_assigned_vertices, candidate_result))
#     return next_vertex_result

# def verify_additional_constraints(self, partially_assigned_vertices: [int], next_vertex: int, potential_mrca: int):
#     """!
#     @brief
#     @param partially_assigned_vertices:
#     @param next_vertex:
#     @param potential_mrca:
#     @return
#     """
#     if len(partially_assigned_vertices) == 2:
#         [first_vertex, second_vertex] = partially_assigned_vertices
#         return self.extend_pair_to_triple(first_vertex, second_vertex, next_vertex, potential_mrca)
#     raise Exception("Unimplemented")

# def extend_pair_to_triple(self, first_vertex: int, second_vertex: int, third_vertex: int,
#                           potential_mrca_candidate: int):
#     """!
#     @brief This function verifies that the Hall condition holds for all the sets involving the third vertex.
#     Before verifying the actual constraint, this function tries to check for the most common scenarios and avoid
#     unnecessary set calculations.
#     @param first_vertex The first vertex from the triplet.
#     @param second_vertex The second vertex from the triplet.
#     @param third_vertex The third vertex from the triplet.
#     @param potential_mrca_candidate A common ancestor of the given triplet of vertices that satisfies the
#     Hall's condition for the (first_vertex, second_vertex) pair.
#     @return True if potential mrca candidate is indeed a potential mrca, False otherwise.
#     """
#     first_vertex_access = self.pedigree.vertex_to_ancestor_map[first_vertex][potential_mrca_candidate]
#     second_vertex_access = self.pedigree.vertex_to_ancestor_map[second_vertex][potential_mrca_candidate]
#     third_vertex_access = self.pedigree.vertex_to_ancestor_map[third_vertex][potential_mrca_candidate]
#     if first_vertex_access == second_vertex_access == third_vertex_access:
#         return False
#     if len(third_vertex_access) == 1:
#         if third_vertex_access == first_vertex_access:
#             return False
#         if third_vertex_access == second_vertex_access:
#             return False
#     return self.triplet_condition_holds(first_vertex, second_vertex, third_vertex, potential_mrca_candidate)
