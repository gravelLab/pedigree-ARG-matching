import os
import random
from pathlib import Path

import networkx
import pytest

from alignment.configuration import MatchingMode, InitialMatchingMode, pedigree_extension
from alignment.graph_matcher import GraphMatcher, get_initial_simulation_mapping_for_mode
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.coalescent_tree import CoalescentTree, SimpleGraph
from itertools import combinations
from collections import Counter

from scripts.utility import get_paths_from_tree_pedigree_directory, dict_has_duplicate_values

logs_folder_name = "logs"
pedigrees_main_folder_name = "pedigrees"


# ----------------------NetworkX implementation--------------------------------------

def verify_valid_alignment_networkx(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
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
    root_vertex_children_number = len(coalescent_tree.children_map[root_vertex])
    assert root_vertex_children_number > 0
    root_vertex_pedigree = alignment[root_vertex]
    network_graph.add_edge(root_vertex_pedigree, target_vertex_label,
                           capacity=root_vertex_children_number)
    proband_number = 0
    for parent in alignment:
        parent_pedigree = alignment[parent]
        if parent not in coalescent_tree.children_map:
            network_graph.add_edge(source_vertex_label, parent_pedigree, capacity=1)
            proband_number += 1
            continue
        children = coalescent_tree.children_map[parent]
        assert len(children) > 0
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
        for parent in graph.parents_map.get(current, []):
            if parent not in paths:
                paths[parent] = []
            paths[parent].append(list(path))
            # Recurse to explore further ancestors
            dfs(parent, path)
        # Backtrack to explore other paths
        path.pop()

    dfs(descendant, [])
    return paths


def build_descendant_ancestor_paths_dict(graph: SimpleGraph) -> dict[int, dict[int, [[int]]]]:
    """
    Build a two-level dictionary where each descendant maps to a dictionary of ancestors
    and their corresponding paths.

    Returns:
        A two-level dictionary:
        - First level key: descendant vertex.
        - Second level key: ancestor vertex.
        - Value: list of all paths from the descendant to the ancestor.
    """
    all_descendant_paths = {vertex: find_paths_for_one_vertex(graph, vertex) for vertex in graph.parents_map}
    return all_descendant_paths


def verify_pedigree_coalescent_tree_alignment(pedigree: PotentialMrcaProcessedGraph, coalescent_tree: CoalescentTree):
    initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=coalescent_tree)
    graph_matcher: GraphMatcher = GraphMatcher(processed_graph=pedigree,
                                               coalescent_tree=coalescent_tree,
                                               matching_mode=MatchingMode.ALL_ALIGNMENTS,
                                               initial_mapping=initial_mapping,
                                               logs_path=None
                                               )
    alignments = graph_matcher.find_mapping()
    coalescent_tree_roots = coalescent_tree.levels[-1]
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
            assert verify_valid_alignment_networkx(
                coalescent_tree=coalescent_tree,
                pedigree=pedigree,
                alignment=alignment,
                root_vertex=coalescent_tree_root
            )


def find_vertex_ancestors(pedigree: PotentialMrcaProcessedGraph, vertex: int):
    vertex_ancestors = set()
    current_level_ancestors = [vertex]
    while current_level_ancestors:
        current_level_ancestors = [vertex for parents in [pedigree.parents_map[x] for x in current_level_ancestors
                                                          if x in pedigree.parents_map]
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
            coalescent_vertex_children = coalescent_tree.children_map.get(coalescent_vertex, [])
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
    for i in range(len(pedigree.levels)):
        for vertex in pedigree.levels[i]:
            assert i == pedigree.vertex_to_level_map[vertex]
            if vertex not in pedigree.children_map:
                continue
            max_child_level = -1
            for child in pedigree.children_map[vertex]:
                assert child in pedigree.vertex_to_level_map
                max_child_level = max(max_child_level, pedigree.vertex_to_level_map[child])
            assert i - max_child_level == 1
    # Test that every non-proband vertex has children
    for i in range(1, len(pedigree.levels)):
        for vertex in pedigree.levels[i]:
            assert vertex in pedigree.children_map


def verify_vertex_to_ancestor_map(pedigree: PotentialMrcaProcessedGraph, vertex_to_ancestors_paths_dict):
    # Test that every non-founder vertex is present in the map and has ancestors
    for level in pedigree.levels[:-1]:
        for vertex in level:
            assert vertex in pedigree.vertex_to_ancestor_map
            ancestors = frozenset(pedigree.vertex_to_ancestor_map[vertex])
            # If the vertex has parents, it must have an ancestor set. Otherwise, the ancestors set must be empty
            if pedigree.parents_map.get(vertex, []):
                ancestors_test = {x for x in vertex_to_ancestors_paths_dict[vertex]}
                assert ancestors == ancestors_test
            else:
                assert not ancestors
            for ancestor in ancestors:
                # Every vertex that is an ancestor to somebody is not a leaf-vertex
                assert pedigree.vertex_to_level_map[ancestor], "A proband vertex is assigned a non-zero level"
                access_vertices = frozenset(pedigree.vertex_to_ancestor_map[vertex][ancestor])
                access_vertices_test = {path[-1] for path in vertex_to_ancestors_paths_dict[vertex][ancestor]}
                assert access_vertices == access_vertices_test, "The access vertices don't match"
                # Verify that every access vertex is an ancestor of vertex
                for access_vertex in access_vertices:
                    assert access_vertex == vertex or access_vertex in ancestors
                    assert ancestor in pedigree.vertex_to_ancestor_map[access_vertex]
                ancestor_children = frozenset(pedigree.children_map[ancestor])
                for non_access_vertex in ancestor_children.difference(access_vertices):
                    assert non_access_vertex not in ancestors
    # Test that every founder vertex has no ancestors
    for vertex in pedigree.levels[-1]:
        assert not pedigree.vertex_to_ancestor_map[vertex]


def verify_preprocessed_pedigree(pedigree: PotentialMrcaProcessedGraph,
                                 vertex_to_ancestors_paths_dict: dict[int, dict[int, [[int]]]]):
    verify_vertex_level_consistency(pedigree)
    verify_vertex_to_ancestor_map(pedigree, vertex_to_ancestors_paths_dict)


def verify_alignment_on_random_small_coalescent_trees(potential_mrca_graph: PotentialMrcaProcessedGraph,
                                                      vertex_to_ancestors_paths_dict: dict[int, dict[int, [[int]]]]):
    pedigree_proband_ploids = potential_mrca_graph.probands
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
                    else:
                        valid_pedigree_mrca.add(common_ancestor)
                        break
        graph = SimpleGraph()
        root = -1
        graph.add_edge(child=first_proband, parent=root)
        graph.add_edge(child=second_proband, parent=root)
        coalescent_tree = CoalescentTree(graph=graph)
        # Flipping the assignments to verify that the symmetry holds
        initial_mapping = {first_proband: [second_proband], second_proband: [first_proband]}
        matcher = GraphMatcher(processed_graph=potential_mrca_graph, coalescent_tree=coalescent_tree,
                               matching_mode=MatchingMode.ALL_ALIGNMENTS,
                               logs_path=None, initial_mapping=initial_mapping)
        result = matcher.find_mapping()
        valid_candidates = {x[root] for x in result[root]}
        assert valid_candidates == valid_pedigree_mrca


def run_verification_for_pedigree_directory(directory_name):
    pedigree_path = None
    coalescent_trees = []
    for file in os.listdir(directory_name):
        file_path = f"{directory_name}/{file}"
        if file.endswith(pedigree_extension):
            # Found the pedigree file
            if pedigree_path is not None:
                raise Exception(f"Multiple pedigrees in {directory_name}")
            pedigree_path = file_path
        else:
            # Found a coalescent tree
            coalescent_trees.append(file_path)
    if pedigree_path is None:
        raise Exception(f"There is no pedigree file in {directory_name}")
    if not coalescent_trees:
        raise Exception(f"There are no coalescent trees in {directory_name}")
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
    root_vertex = tree.get_root_vertex()
    initial_mapping = get_initial_simulation_mapping_for_mode(coalescent_tree=tree,
                                                              mode=InitialMatchingMode.INDIVIDUAL)
    matcher = GraphMatcher(
        processed_graph=pedigree,
        coalescent_tree=tree,
        matching_mode=MatchingMode.ALL_ALIGNMENTS,
        logs_path=None,
        initial_mapping=initial_mapping
    )
    result = matcher.find_mapping()
    expected_alignments = [
        {2: 2, 4: 4, 3: 6},
        {2: 2, 4: 4, 3: 7},
        {2: 2, 4: 4, 3: 12},
        {2: 2, 4: 4, 3: 13}
    ]
    assert root_vertex in result, "The root vertex hasn't been mapped"
    alignments = result[root_vertex]
    assert dict_list_equal_ignore_order(alignments, expected_alignments), "Didn't find the expected alignments"
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
        assert not matcher.verify_valid_alignment(invalid_alignment)


def find_alignments_for_path(directory_path: str | Path, initial_mapping: dict[int, [int]]):
    paths = get_paths_from_tree_pedigree_directory(directory_path)
    if not paths:
        raise ValueError(f"{directory_path} is not a valid tree-pedigree directory")
    pedigree_path, tree_path = paths
    pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path)
    tree = CoalescentTree.get_coalescent_tree_from_file(tree_path)
    tree_root = tree.get_root_vertex()
    matcher = GraphMatcher(
        processed_graph=pedigree,
        coalescent_tree=tree,
        initial_mapping=initial_mapping,
        matching_mode=MatchingMode.ALL_ALIGNMENTS,
        logs_path=None,
    )
    alignments_map = matcher.find_mapping()
    alignments = alignments_map[tree_root]
    return matcher, alignments


@pytest.fixture
def test_collision_is_rejected_1_path():
    return Path(os.path.dirname(os.path.abspath(__file__))) / "test_collision_is_rejected_1"


def test_collision_is_rejected_1(test_collision_is_rejected_1_path):
    initial_mapping = {
        7: [2],
        8: [4],
        10: [14]
    }
    matcher, alignments = find_alignments_for_path(directory_path=test_collision_is_rejected_1_path,
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
    assert not matcher.verify_valid_alignment(alignment=invalid_alignment)
    invalid_alignment[9] = 13
    assert not matcher.verify_valid_alignment(alignment=invalid_alignment)


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
    _, alignments = find_alignments_for_path(directory_path=test_collision_is_rejected_2_path,
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
