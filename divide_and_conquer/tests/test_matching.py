import os
import unittest

from divide_and_conquer.graph_matcher import GraphMather, MatcherLogger
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from divide_and_conquer.subtree_matcher import SubtreeMatcher
from genealogical_graph import CoalescentTree
from graph import Graph

pedigrees_main_folder_name = "pedigrees"
pedigree_extension = ".pedigree"


# TODO: Consider using subtests

def test_pedigree_coalescent_tree_alignment(pedigree: PotentialMrcaProcessedGraph,
                                            coalescent_tree: CoalescentTree):
    logger = MatcherLogger()
    graph_matcher: GraphMather = GraphMather(pedigree, coalescent_tree, logger)
    alignments = graph_matcher.find_mapping()
    test_identity_is_present(coalescent_tree=coalescent_tree, alignments=alignments)
    test_all_found_alignments_are_correct(pedigree=pedigree, coalescent_tree=coalescent_tree, alignments=alignments)


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


def test_all_found_alignments_are_correct(pedigree: PotentialMrcaProcessedGraph,
                                          coalescent_tree: CoalescentTree,
                                          alignments: {int: [SubtreeMatcher]}):
    coalescent_tree_roots = coalescent_tree.levels[-1]
    for coalescent_tree_root in coalescent_tree_roots:
        test_alignments_correct_for_clade(coalescent_tree=coalescent_tree, pedigree=pedigree,
                                          coalescent_tree_root=coalescent_tree_root, alignments=alignments)


def test_alignments_correct_for_clade(coalescent_tree: CoalescentTree, pedigree: PotentialMrcaProcessedGraph,
                                      coalescent_tree_root: int, alignments: {int: [SubtreeMatcher]}):
    if coalescent_tree_root not in coalescent_tree.children_map:
        return
    coalescent_tree_root_children = coalescent_tree.children_map[coalescent_tree_root]
    # Need to add an if statement because we are using a workaround for this case
    if len(coalescent_tree_root_children) > 1:
        for alignment in alignments[coalescent_tree_root]:
            alignment: SubtreeMatcher
            if alignment.children_assignments is None:
                return
            children_vertices = list(alignment.children_assignments.values())
            children_vertices_pmracs = find_pmrcas_for_vertices(pedigree, children_vertices)
            if alignment.root_pedigree not in children_vertices_pmracs:
                raise Exception("Assigned coalescent vertex is not a MRCA for the children assignments")
    for coalescent_tree_root_child in coalescent_tree_root_children:
        test_alignments_correct_for_clade(coalescent_tree=coalescent_tree, pedigree=pedigree,
                                          coalescent_tree_root=coalescent_tree_root_child, alignments=alignments)


def test_identity_is_present(coalescent_tree: CoalescentTree, alignments: {int: [SubtreeMatcher]}):
    coalescent_tree_roots = coalescent_tree.levels[-1]
    for coalescent_tree_root in coalescent_tree_roots:
        test_identity_present_for_clade(coalescent_tree=coalescent_tree,
                                        coalescent_tree_root=coalescent_tree_root,
                                        alignments=alignments)


def test_identity_present_for_clade(coalescent_tree: CoalescentTree,
                                    coalescent_tree_root: int, alignments: {int: [SubtreeMatcher]}):
    root_identity_alignments = [x for x in alignments[coalescent_tree_root]
                                if x.root_pedigree == x.root_coalescent_tree]
    if not root_identity_alignments:
        raise Exception("No identity alignments")
    # If coalescent_tree_root is a proband in the coalescent tree, there is nothing else to verify
    if coalescent_tree_root not in coalescent_tree.children_map:
        return
    # Otherwise, check that the children subtrees also have an identity mapping
    identity_present = False
    children_identity_alignment = {x: x for x in coalescent_tree.children_map[coalescent_tree_root]}
    correct_children_set = set(coalescent_tree.children_map[coalescent_tree_root])
    for alignment in root_identity_alignments:
        alignment: SubtreeMatcher
        children_set = set(alignment.children_assignments.values())
        if children_set == correct_children_set:
            identity_present = True
            break
    if not identity_present:
        raise Exception("There are no identity mappings for the children")
    for child in coalescent_tree.children_map[coalescent_tree_root]:
        test_identity_present_for_clade(coalescent_tree=coalescent_tree,
                                        coalescent_tree_root=child,
                                        alignments=alignments)


def test_vertex_level_consistency(pedigree: PotentialMrcaProcessedGraph):
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


def test_vertex_to_ancestor_map(pedigree: PotentialMrcaProcessedGraph):
    levels_number = len(pedigree.levels)
    # Test that every non-founder vertex is present in the map and has ancestors
    for i in range(levels_number - 1):
        for vertex in pedigree.levels[i]:
            assert vertex in pedigree.vertex_to_ancestor_map
            ancestors = set(pedigree.vertex_to_ancestor_map[vertex])
            ancestors_test = find_vertex_ancestors(pedigree, vertex)
            assert ancestors == ancestors_test
            for ancestor in ancestors:
                assert pedigree.vertex_to_level_map[ancestor]
                ancestor_children = set(pedigree.children_map[ancestor])
                access_vertices = set(pedigree.vertex_to_ancestor_map[vertex][ancestor])
                # Verify that the access vertices are the ancestor children
                assert access_vertices <= ancestor_children
                assert access_vertices
                # Verify that every access vertex is an ancestor of vertex
                for access_vertex in access_vertices:
                    assert access_vertex == vertex or access_vertex in ancestors
                for non_access_vertex in ancestor_children.difference(access_vertices):
                    assert non_access_vertex not in ancestors
    # Test that every founder vertex has no ancestors
    for vertex in pedigree.levels[levels_number - 1]:
        assert not len(pedigree.vertex_to_ancestor_map[vertex])


def test_preprocessed_pedigree(pedigree: PotentialMrcaProcessedGraph):
    test_vertex_level_consistency(pedigree)
    test_vertex_to_ancestor_map(pedigree)


def test_pedigree_directory(directory_name):
    pedigree = None
    coalescent_trees = []
    for file in os.listdir(directory_name):
        file_path = f"{directory_name}/{file}"
        if file.endswith(pedigree_extension):
            # Found the pedigree file
            if pedigree is not None:
                raise Exception(f"Multiple pedigrees in {directory_name}")
            pedigree = file_path
        else:
            # Found a coalescent tree
            coalescent_trees.append(file_path)
    if pedigree is None:
        raise Exception(f"There is no pedigree file in {directory_name}")
    if not coalescent_trees:
        raise Exception(f"There are no coalescent trees in {directory_name}")
    potential_mrca_graph = PotentialMrcaProcessedGraph(pedigree=Graph.get_pedigree_from_file(filename=pedigree))
    test_preprocessed_pedigree(potential_mrca_graph)
    for coalescent_tree_name in coalescent_trees:
        # Verify that all the alignments are consistent
        coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(coalescent_tree_name)
        test_pedigree_coalescent_tree_alignment(pedigree=potential_mrca_graph, coalescent_tree=coalescent_tree)
        print(f"Alignment with {coalescent_tree_name} has been tested")


class TestPotentialCommonAncestors(unittest.TestCase):

    def test_generated_pedigrees_alignments(self):
        os.chdir(pedigrees_main_folder_name)
        for directory_entry in os.listdir():
            if os.path.isdir(directory_entry):
                print(f"Testing the algorithm on {directory_entry}")
                test_pedigree_directory(directory_name=directory_entry)
                print(f"{directory_entry} has been tested")


if __name__ == '__main__':
    unittest.main()
