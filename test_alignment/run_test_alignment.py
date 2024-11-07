import os

from alignment.graph_matcher import MatcherLogger, GraphMatcher
from alignment.potential_mrca_processed_graph import PotentialMrcaProcessedGraph
from graph.coalescent_tree import CoalescentTree
from scripts.run_alignment import save_alignment_result_to_files
from scripts.utility import get_file_path

# initial_mapping = {
#     0: [99539940, 99539941],
#     1: [99539940, 99539941],
#     2: [99575044, 99575045],
#     3: [99575044, 99575045],
#     4: [99955730, 99955731],
#     5: [99955730, 99955731],
#     6: [100077508, 100077509],
#     7: [100077508, 100077509],
#     8: [100170062, 100170063],
#     9: [100170062, 100170063],
#     10: [100172634, 100172635],
#     11: [100172634, 100172635],
# }

initial_mapping = {
    0: [99539940, 99539941],
    2: [99575044, 99575045],
    4: [99955730, 99955731],
    7: [100077508, 100077509],
    8: [100170062, 100170063],
    11: [100172634, 100172635],
}

pedigree_path = get_file_path("Specify the path to the pedigree:")
tree_path = get_file_path("Specify the path to the tree:")
coalescent_tree_filename = os.path.basename(tree_path)
coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path)
coalescent_tree.reduce_to_ascending_genealogy(probands=initial_mapping.keys())
coalescent_tree.remove_unary_nodes()
pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path, separation_symbol='\t',
                                                                     missing_parent_notation=["0"],
                                                                     initialize_levels=False,
                                                                     initialize_ancestor_maps=False
                                                                     )
pedigree_probands = {y for x in initial_mapping.values() for y in x}
pedigree.reduce_to_ascending_genealogy(probands=pedigree_probands, recalculate_levels=True)
pedigree.initialize_potential_mrca_map()
logger = MatcherLogger()
graph_matcher = GraphMatcher(processed_graph=pedigree, coalescent_tree=coalescent_tree,
                             logger=logger, initial_mapping=initial_mapping)
results = graph_matcher.find_mapping()
os.mkdir("results")
os.chdir("results")
save_alignment_result_to_files(alignment_result=results, coalescent_tree=coalescent_tree,
                               pedigree=pedigree, coalescent_tree_filename=coalescent_tree_filename)
