import os
import random
from pathlib import Path
from typing import cast

import networkx
import numpy
import pytest
from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.GenGraph import GenGraph

from alignment.alignment_result import FullAlignmentResult, AlignmentResult, FailedClimbingAlignmentResult, \
    SuccessCladeAlignmentResults
from alignment.configuration import AlignmentVertexMode, ProbandInitialAssignmentsMode, AlignmentEdgeMode
from alignment.graph_matcher import GraphMatcher, get_initial_simulation_mapping_for_mode
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from itertools import combinations
from collections import Counter

from scripts.utility.basic_utility import get_paths_from_tree_pedigree_directory
from scripts.utility.alignment_utility import dict_has_duplicate_values

logs_folder_name = "logs"
pedigrees_main_folder_name = "pedigrees"


# ----------------------NetworkX implementation--------------------------------------

def verify_alignment_necessary_condition_networkx(coalescent_tree: CoalescentTree,
                                                  pedigree: PotentialMrcaProcessedGraph,
                                                  alignment: dict[int, int], root_vertex: int):
    """
    Verifies that the given alignment is correct.

    Args:
        coalescent_tree: The coalescent tree
        pedigree: The pedigree
        alignment (dict): The alignment to be verified.
        root_vertex (int): The root vertex of the aligned clade.

    Returns:
        Whether the given alignment is correct.
    """
    if dict_has_duplicate_values(alignment):
        return False
    network_graph = networkx.DiGraph()
    source_vertex_label = "s"
    target_vertex_label = "t"
    # Adding the edge from the root to the sink vertex
    root_vertex_children_number = len(coalescent_tree.get_children(root_vertex))
    if not root_vertex_children_number:
        return alignment[root_vertex] in pedigree
    root_vertex_pedigree = alignment[root_vertex]
    network_graph.add_edge(root_vertex_pedigree, target_vertex_label,
                           capacity=root_vertex_children_number)
    proband_number = 0
    for parent in alignment:
        parent_pedigree = alignment[parent]
        children = coalescent_tree.get_children(parent)
        if not children:
            network_graph.add_edge(source_vertex_label, parent_pedigree, capacity=1)
            proband_number += 1
            continue
        children_pedigree = [alignment[x] for x in children if alignment.get(x) is not None]
        children_number = len(children_pedigree)
        add_edges_to_mrca_from_descendants_networkx(pedigree, network_graph, parent_pedigree, children_pedigree)
        if parent_pedigree != root_vertex_pedigree:
            network_graph.add_edge(parent_pedigree, target_vertex_label, capacity=children_number - 1)
    maximum_flow = networkx.maximum_flow_value(flowG=network_graph, _s=source_vertex_label, _t=target_vertex_label)
    return maximum_flow == proband_number


def add_edges_to_mrca_from_descendants_networkx(pedigree: PotentialMrcaProcessedGraph,
                                                flow_network: networkx.DiGraph, mrca: int,
                                                descendant_vertices: [int]):
    """
    Adds all the paths from the given descendant vertices to the specified mrca to the specified graph.

    Args:
        pedigree: The pedigree.
        flow_network (networkx.DiGraph): The graph to which the edges are added.
        mrca (int): The mrca of the descendant vertices.
        descendant_vertices (list[int]): The descendant vertices sharing the mrca.
    """
    mrca_access_label = f"{mrca}'"
    flow_network.add_edge(mrca_access_label, mrca, capacity=len(descendant_vertices))
    for descendant in descendant_vertices:
        current_level_vertices = [mrca]
        while current_level_vertices:
            next_level_vertices = set()
            for vertex in current_level_vertices:
                if vertex == descendant:
                    continue
                neighbors = pedigree.vertex_to_ancestor_map[descendant].get(vertex, [])
                for neighbor in neighbors:
                    vertex_access_label = f"{vertex}'"
                    # We can potentially override the edge's capacity if vertex is the mrca in the
                    # coalescent tree
                    if not flow_network.has_edge(vertex_access_label, vertex):
                        flow_network.add_edge(vertex_access_label, vertex, capacity=1)
                    flow_network.add_edge(neighbor, vertex_access_label, capacity=1)
                next_level_vertices.update(neighbors)
            current_level_vertices = next_level_vertices


# --------------------------- End NetworkX implementation-----------------------------------

def dict_list_equal_ignore_order(list1, list2):
    return Counter(frozenset(d.items()) for d in list1) == Counter(frozenset(d.items()) for d in list2)


def find_paths_for_one_vertex(graph: PotentialMrcaProcessedGraph, descendant: int):
    """
    Find all paths from the given descendant to its ancestors.
    """
    paths = {}

    def dfs(current, path):
        # Add the current node to the path
        path.append(current)
        for parent in graph.get_parents(current):
            if parent not in paths:
                paths[parent] = []
            paths[parent].append(list(path))
            # Recurse to explore further ancestors
            dfs(parent, path)
        # Backtrack to explore other paths
        path.pop()

    dfs(descendant, [])
    return paths


def build_descendant_ancestor_paths_dict(graph: GenGraph) -> dict[int, dict[int, [[int]]]]:
    """
    Build a two-level dictionary where each descendant maps to a dictionary of ancestors
    and their corresponding paths.

    Returns:
        A two-level dictionary:
        - First level key: descendant vertex.
        - Second level key: ancestor vertex.
        - Value: list of all paths from the descendant to the ancestor.
    """
    all_descendant_paths = {vertex: find_paths_for_one_vertex(graph, vertex) for vertex in graph}
    return all_descendant_paths


def verify_pedigree_coalescent_tree_alignment(pedigree: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree):
    initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
    alignments = []

    def collect_results(alignment_result: FullAlignmentResult):
        assert alignment_result.is_valid
        alignments.append(alignment_result.vertex_alignment)

    graph_matcher: GraphMatcher = GraphMatcher(processed_graph=pedigree,
                                               coalescent_tree=coalescent_tree,
                                               alignment_vertex_mode=AlignmentVertexMode.ALL_ALIGNMENTS,
                                               alignment_edge_mode=AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT,
                                               initial_mapping=initial_mapping,
                                               logs_path=None,
                                               result_callback_function=collect_results,
                                               calculate_posterior_probabilities=False
                                               )
    graph_matcher.find_alignments()
    coalescent_tree_roots = coalescent_tree.get_founders()
    for coalescent_tree_root in coalescent_tree_roots:
        clade_alignments = alignments[coalescent_tree_root]
        assert clade_alignments
        verify_identity_present_for_clade(coalescent_tree=coalescent_tree,
                                          coalescent_tree_root=coalescent_tree_root,
                                          clade_alignments=clade_alignments)
        verify_clade_alignments_locally_consistent(coalescent_tree=coalescent_tree, pedigree=pedigree,
                                                   clade_alignments=clade_alignments)
        verify_all_root_candidate_ploids_found(clade_alignments=clade_alignments, clade_root=coalescent_tree_root)
        for alignment in clade_alignments:
            assert verify_alignment_necessary_condition_networkx(
                coalescent_tree=coalescent_tree,
                pedigree=pedigree,
                alignment=alignment,
                root_vertex=coalescent_tree_root
            )


def find_vertex_ancestors(pedigree: PotentialMrcaProcessedGraph, vertex: int):
    vertex_ancestors = set()
    current_level_ancestors = [vertex]
    while current_level_ancestors:
        current_level_ancestors = [vertex for parents in [pedigree.get_parents(x) for x in current_level_ancestors]
                                   for vertex in parents]
        vertex_ancestors.update(current_level_ancestors)
    return vertex_ancestors


# TODO: Add caching to speed up the tests (read the framework documentation regarding subtests, setups and clean-ups)
def find_pmrcas_for_vertices(pedigree: PotentialMrcaProcessedGraph, vertices: [int]):
    vertices_ancestors_list = list()
    for vertex in vertices:
        vertices_ancestors_list.append(find_vertex_ancestors(pedigree, vertex))
    return set.intersection(*vertices_ancestors_list)


def halls_condition_holds(mrca: int, descendants: [int], vertex_to_ancestor_map: dict[int, dict[int, tuple[int]]]):
    """
    Checks if Hall's condition holds for the given set of vertices.
    """
    n = len(descendants)
    # Generate all possible subsets of the vertex set (excluding the empty set)
    for r in range(1, n + 1):
        for subset in combinations(descendants, r):
            neighbor_set = set()
            for vertex in subset:
                neighbor_set.update(vertex_to_ancestor_map[vertex][mrca])
            if len(neighbor_set) < len(subset):
                return False
    return True


def verify_clade_alignments_locally_consistent(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                                               clade_alignments: [dict[int, int]]):
    for alignment in clade_alignments:
        for coalescent_vertex, vertex_pedigree_assignment in alignment.items():
            coalescent_vertex_children = coalescent_tree.get_children(coalescent_vertex)
            if not coalescent_vertex_children:
                # Picked a proband vertex, there is nothing to be verified, so we can simply continue
                continue
            children_pedigree_assignments = [alignment[x] for x in coalescent_vertex_children]
            for children_assignment in children_pedigree_assignments:
                assert vertex_pedigree_assignment in pedigree.vertex_to_ancestor_map[children_assignment]
            assert halls_condition_holds(vertex_pedigree_assignment, children_pedigree_assignments,
                                         pedigree.vertex_to_ancestor_map), \
                (f"The inferred coalescent tree assignment {vertex_pedigree_assignment} is not"
                 f" a pmrca of {children_pedigree_assignments} in the alignment {alignment}")


def verify_all_root_candidate_ploids_found(clade_alignments: [dict[int, int]], clade_root: int):
    # The alignment problem has a natural symmetry. Specifically, if we infer that a specific pedigree
    # ploid is a valid assignment for the clade's root, then the other ploid of the same individual
    # is also a valid assignment! Indeed, the other ploid is connected to the same ploids as the initial ploid.
    # This test verifies that this holds for all the found alignments
    pedigree_ploid_candidates = {clade_alignment[clade_root] for clade_alignment in clade_alignments}
    pedigree_individual_candidates = {x // 2 for x in pedigree_ploid_candidates}
    respective_ploid_candidates = {val for x in pedigree_individual_candidates for val in (2 * x, 2 * x + 1)}
    assert pedigree_ploid_candidates == respective_ploid_candidates


def verify_identity_present_for_clade(coalescent_tree: CoalescentTree, coalescent_tree_root: int,
                                      clade_alignments: [dict[int, int]]):
    clade_vertices = coalescent_tree.get_connected_component_for_vertex(coalescent_tree_root)
    identity_solution = {x: x for x in clade_vertices}
    assert identity_solution in clade_alignments, "The identity solution is not present in the solution set"


def verify_vertex_level_consistency(pedigree: PotentialMrcaProcessedGraph):
    # Test that the level assignments are consistent
    pedigree_levels = pedigree.get_levels()
    for i, level in enumerate(pedigree_levels):
        for vertex in level:
            assert i == pedigree.get_vertex_level(vertex)
            vertex_children = pedigree.get_children(vertex)
            if not vertex_children:
                continue
            max_child_level = -1
            for child in vertex_children:
                max_child_level = max(max_child_level, pedigree.get_vertex_level(child))
            assert i - max_child_level == 1
    # Test that every non-proband vertex has children
    for level in pedigree_levels[1:]:
        for vertex in level:
            assert pedigree.has_children(vertex)


def verify_vertex_to_ancestor_map(pedigree: PotentialMrcaProcessedGraph, vertex_to_ancestors_paths_dict):
    # Test that every non-founder vertex is present in the map and has ancestors
    pedigree_levels = pedigree.get_levels()
    for level in pedigree_levels[:-1]:
        for vertex in level:
            assert vertex in pedigree.vertex_to_ancestor_map
            ancestors = frozenset(pedigree.vertex_to_ancestor_map[vertex])
            # If the vertex has parents, it must have an ancestor set. Otherwise, the ancestors set must be empty
            if pedigree.has_parents(vertex):
                ancestors_test = {x for x in vertex_to_ancestors_paths_dict[vertex]}
                assert ancestors == ancestors_test
            else:
                assert not ancestors
            for ancestor in ancestors:
                # Every vertex that is an ancestor to somebody is not a leaf-vertex
                assert pedigree.get_vertex_level(ancestor), "A proband vertex is assigned a non-zero level"
                access_vertices = frozenset(pedigree.vertex_to_ancestor_map[vertex][ancestor])
                access_vertices_test = {path[-1] for path in vertex_to_ancestors_paths_dict[vertex][ancestor]}
                assert access_vertices == access_vertices_test, "The access vertices don't match"
                # Verify that every access vertex is an ancestor of vertex
                for access_vertex in access_vertices:
                    assert access_vertex == vertex or access_vertex in ancestors
                    assert ancestor in pedigree.vertex_to_ancestor_map[access_vertex]
                ancestor_children = frozenset(pedigree.get_children(ancestor))
                for non_access_vertex in ancestor_children.difference(access_vertices):
                    assert non_access_vertex not in ancestors
    # Test that every founder vertex has no ancestors
    for vertex in pedigree_levels[-1]:
        assert not pedigree.vertex_to_ancestor_map[vertex]


def verify_preprocessed_pedigree(pedigree: PotentialMrcaProcessedGraph,
                                 vertex_to_ancestors_paths_dict: dict[int, dict[int, [[int]]]]):
    verify_vertex_level_consistency(pedigree)
    verify_vertex_to_ancestor_map(pedigree, vertex_to_ancestors_paths_dict)


def verify_alignment_on_random_small_coalescent_trees(potential_mrca_graph: PotentialMrcaProcessedGraph,
                                                      vertex_to_ancestors_paths_dict: dict[int, dict[int, [[int]]]]):
    pedigree_proband_ploids = potential_mrca_graph.get_sink_vertices()
    for i in range(10):
        tree_probands = list(random.sample(sorted(pedigree_proband_ploids), 2))
        [first_proband, second_proband] = tree_probands
        # Finding valid pedigree MRCA for these vertices
        common_ancestors = set.intersection(*[set(vertex_to_ancestors_paths_dict[x]) for x in tree_probands])
        valid_pedigree_mrca = set()
        for common_ancestor in common_ancestors:
            vertex_to_paths = {x: vertex_to_ancestors_paths_dict[x][common_ancestor] for x in tree_probands}
            # Check if the paths intersect only at the common ancestor
            for first_path in vertex_to_paths[tree_probands[0]]:
                for second_path in vertex_to_paths[tree_probands[1]]:
                    path1_set = frozenset(first_path)
                    path2_set = frozenset(second_path)
                    # If the intersection is not empty, paths intersect before the common ancestor
                    if path1_set & path2_set:
                        continue
                    valid_pedigree_mrca.add(common_ancestor)
                    break
        coalescent_tree = CoalescentTree()
        root = -1
        coalescent_tree.add_edge(child=first_proband, parent=root)
        coalescent_tree.add_edge(child=second_proband, parent=root)
        # Flipping the assignments to verify that the symmetry holds
        initial_mapping = {first_proband: [second_proband], second_proband: [first_proband]}
        result = []

        def collect_results(alignment_result: AlignmentResult):
            match alignment_result:
                case FailedClimbingAlignmentResult():
                    return
                case FullAlignmentResult():
                    alignment_result = cast(FullAlignmentResult, alignment_result)
                    assert alignment_result.is_valid
                    assert alignment_result.clade_root is not None
                    if alignment_result.clade_root == root:
                        result.append(alignment_result.vertex_alignment)

        matcher = GraphMatcher(processed_graph=potential_mrca_graph, coalescent_tree=coalescent_tree,
                               alignment_vertex_mode=AlignmentVertexMode.ALL_ALIGNMENTS,
                               alignment_edge_mode=AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT,
                               logs_path=None, initial_mapping=initial_mapping,
                               result_callback_function=collect_results,
                               calculate_posterior_probabilities=False)
        matcher.find_alignments()
        valid_candidates = {x[root] for x in result}
        assert valid_candidates == valid_pedigree_mrca


def run_verification_for_pedigree_directory(directory_name):
    coalescent_trees = []
    paths = get_paths_from_tree_pedigree_directory(directory_name)
    if not paths:
        raise ValueError(f"{invalid_alignment_discarded_1} is not a valid tree-pedigree directory")
    pedigree_path, tree_path = paths
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path)
    vertex_to_ancestors_paths_dict = build_descendant_ancestor_paths_dict(pedigree)
    verify_preprocessed_pedigree(pedigree, vertex_to_ancestors_paths_dict)
    for coalescent_tree_name in coalescent_trees:
        # Verify that all the alignments are "locally" consistent. Specifically, we want to verify that
        # every non-leaf vertex in the coalescent tree is assigned to a pedigree ploid that is a valid pmrca
        # for the children vertices
        coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(coalescent_tree_name)
        coalescent_tree.remove_unary_nodes()
        verify_pedigree_coalescent_tree_alignment(pedigree=pedigree, coalescent_tree=coalescent_tree)
        print(f"Alignment with {coalescent_tree_name} has been tested")
    verify_alignment_on_random_small_coalescent_trees(pedigree, vertex_to_ancestors_paths_dict)


def get_pedigree_directories():
    os.chdir(pedigrees_main_folder_name)
    return [
        directory_entry
        for directory_entry in os.listdir()
        if os.path.isdir(directory_entry) and directory_entry != logs_folder_name
    ]


@pytest.mark.parametrize("directory_entry", get_pedigree_directories())
def test_generated_pedigrees_alignments(directory_entry):
    print(f"Testing the algorithm on {directory_entry}")
    run_verification_for_pedigree_directory(directory_name=directory_entry)
    print(f"{directory_entry} has been tested")


@pytest.fixture
def invalid_alignment_discarded_1():
    return Path(os.path.dirname(os.path.abspath(__file__))) / "test_invalid_alignment_is_discarded_1"


def test_invalid_alignment_discarded_1(invalid_alignment_discarded_1):
    paths = get_paths_from_tree_pedigree_directory(invalid_alignment_discarded_1)
    if not paths:
        raise ValueError(f"{invalid_alignment_discarded_1} is not a valid tree-pedigree directory")
    pedigree_path, tree_path = paths
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(pedigree_path)
    tree = CoalescentTree.get_coalescent_tree_from_file(tree_path)
    initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=tree,
                                                              mode=ProbandInitialAssignmentsMode.INDIVIDUAL)
    alignments: [FullAlignmentResult] = []

    def collect_results(alignment_result: FullAlignmentResult):
        assert alignment_result.is_valid
        alignments.append(alignment_result)

    matcher = GraphMatcher(
        processed_graph=pedigree,
        coalescent_tree=tree,
        alignment_vertex_mode=AlignmentVertexMode.ALL_ALIGNMENTS,
        alignment_edge_mode=AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS,
        logs_path=None,
        initial_mapping=initial_mapping,
        result_callback_function=collect_results,
        calculate_posterior_probabilities=True
    )
    matcher.find_alignments()
    clade_results: SuccessCladeAlignmentResults = SuccessCladeAlignmentResults(
        clade_root=tree.get_root_vertex(),
        alignments=alignments
    )
    clade_results.calculate_posterior_probabilities()
    alignment_to_posterior_probability = clade_results.clade_alignment_posterior_probabilities.vertex_alignment_to_posterior_probability
    vertex_to_posterior_probability = clade_results.clade_alignment_posterior_probabilities.vertex_posterior_probabilities
    expected_alignments = [
        {2: 2, 4: 4, 3: 6},
        {2: 2, 4: 4, 3: 7},
        {2: 2, 4: 4, 3: 12},
        {2: 2, 4: 4, 3: 13}
    ]
    alignment_posterior_probabilities = [4/9, 4/9, 1/18, 1/18]
    alignment_simple_probabilities = [0.25, 0.25, 2 ** -5, 2 ** -5]
    expected_vertex_to_posterior_probability = {
        2: 1, 4: 1, 6: 5/9, 7: 5/9, 12: 1/18, 13: 1/18,
        8: 1/9, 10: 1/9,
    }
    resulting_vertex_alignments = [x.vertex_alignment for x in alignments]
    assert dict_list_equal_ignore_order(resulting_vertex_alignments, expected_alignments), \
        "Didn't find the expected alignments"
    invalid_alignments = [
        {
            2: 6, 4: 8, 3: 12
        },
        {
            2: 6, 4: 8, 3: 13
        },
        {
            2: 7, 4: 10, 3: 12
        },
        {
            2: 7, 4: 10, 3: 13
        },
        {
            2: 9, 4: 10, 3: 13
        },
        {
            2: 7, 4: 11, 3: 13
        },
        {
            2: 7, 4: 8, 3: 2
        }
    ]
    for invalid_alignment in invalid_alignments:
        assert not matcher._verify_valid_alignment(invalid_alignment)
    # Verify the posterior probabilities
    for alignment in alignments:
        alignment: FullAlignmentResult
        assert alignment.posterior_probabilities
        vertex_alignment = alignment.vertex_alignment
        try:
            index = expected_alignments.index(vertex_alignment)
        except ValueError:
            pytest.fail(f"Received an unexpected alignment: {vertex_alignment}")
        expected_alignment_probability = alignment_simple_probabilities[index]
        actual_alignment_probability = float(alignment.posterior_probabilities.vertex_alignment_probability)
        assert numpy.isclose(actual_alignment_probability, expected_alignment_probability), \
            (f"Invalid vertex alignment probability {actual_alignment_probability}, "
             f"expected {expected_alignment_probability} instead. The alignment: {vertex_alignment}")
        expected_alignment_posterior_probability = alignment_posterior_probabilities[index]
        actual_alignment_posterior_probability = alignment_to_posterior_probability[index]
        assert numpy.isclose(actual_alignment_posterior_probability, expected_alignment_posterior_probability), \
            (f"Invalid vertex alignment posterior probability {actual_alignment_probability}, "
             f"expected {expected_alignment_probability} instead. The alignment: {vertex_alignment}")
    assert len(vertex_to_posterior_probability) == len(expected_vertex_to_posterior_probability)
    for vertex, vertex_posterior_probability in expected_vertex_to_posterior_probability.items():
        result_posterior_probability = vertex_to_posterior_probability[vertex]
        assert numpy.isclose(vertex_posterior_probability, result_posterior_probability), \
            (f"Invalid posterior probability for vertex {vertex}: {result_posterior_probability}, "
             f"expected {vertex_posterior_probability} instead")


def find_alignments_for_tree_pedigree_directory(directory_path: str | Path, initial_mapping: dict[int, [int]]):
    paths = get_paths_from_tree_pedigree_directory(directory_path)
    if not paths:
        raise ValueError(f"{directory_path} is not a valid tree-pedigree directory")
    pedigree_path, tree_path = paths
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path)
    tree = CoalescentTree.get_coalescent_tree_from_file(tree_path)
    result = []

    def collect_results(alignment_result: FullAlignmentResult):
        assert alignment_result.is_valid
        # if tree_root in alignment_result.vertex_alignment:
        result.append(alignment_result.vertex_alignment)

    matcher = GraphMatcher(
        processed_graph=pedigree,
        coalescent_tree=tree,
        initial_mapping=initial_mapping,
        alignment_vertex_mode=AlignmentVertexMode.ALL_ALIGNMENTS,
        alignment_edge_mode=AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT,
        logs_path=None,
        result_callback_function=collect_results,
        calculate_posterior_probabilities=False
    )
    matcher.find_alignments()
    return matcher, result


@pytest.fixture
def test_collision_is_rejected_1_path():
    return Path(os.path.dirname(os.path.abspath(__file__))) / "test_collision_is_rejected_1"


def test_collision_is_rejected_1(test_collision_is_rejected_1_path):
    initial_mapping = {
        7: [2],
        8: [4],
        10: [14]
    }
    matcher, alignments = find_alignments_for_tree_pedigree_directory(directory_path=test_collision_is_rejected_1_path,
                                                                      initial_mapping=initial_mapping
                                                                      )
    expected_alignments = [
        {
            7: 2,
            8: 4,
            10: 14,
            9: 6
        },
        {
            7: 2,
            8: 4,
            10: 14,
            9: 7
        }
    ]
    assert dict_list_equal_ignore_order(expected_alignments, alignments)
    invalid_alignment = {
        7: 2,
        8: 4,
        10: 14,
        9: 12
    }
    assert not matcher._verify_valid_alignment(potential_alignment=invalid_alignment)
    invalid_alignment[9] = 13
    assert not matcher._verify_valid_alignment(potential_alignment=invalid_alignment)


@pytest.fixture
def test_collision_is_rejected_2_path():
    return Path(os.path.dirname(os.path.abspath(__file__))) / "test_collision_is_rejected_2"


def test_collision_is_rejected_2(test_collision_is_rejected_2_path):
    initial_mapping = {
        1: [2],
        2: [4],
        4: [14],
        5: [16],
    }
    _, alignments = find_alignments_for_tree_pedigree_directory(directory_path=test_collision_is_rejected_2_path,
                                                                initial_mapping=initial_mapping)
    expected_alignments = [
        {
            1: 2,
            2: 4,
            4: 14,
            5: 16,
            3: 6,
            6: 18,
            7: 20
        },
        {
            1: 2,
            2: 4,
            4: 14,
            5: 16,
            3: 6,
            6: 18,
            7: 21
        }
    ]
    assert dict_list_equal_ignore_order(expected_alignments, alignments)
