import os
import time

from genealogical_graph import *
from divide_and_conquer.graph_matcher import GraphMather, MatcherLogger
from divide_and_conquer.potential_mrca_processed_graph import PotentialMrcaProcessedGraph

directory = input("Specify the testing directory in the pedigrees folder:")
os.chdir(f"pedigrees/{directory}")

# filename_input = input("Specify the pedigree file:")
filename = f"{directory}.pedigree"

pedigree = Graph.get_pedigree_from_file(filename=filename)
start_preprocessing = time.time()
potential_mrca_graph = PotentialMrcaProcessedGraph(pedigree=pedigree)
end_preprocessing = time.time()

print(f"Preprocessing time: {end_preprocessing - start_preprocessing} seconds")

# Get the current working directory
current_directory = os.getcwd()
# Get a list of all files in the current directory
total_alignment_time = 0
count = 0
for file in os.listdir(current_directory):
    filename = os.path.basename(file)
    if filename == f"{directory}.pedigree":
        continue
    print(f"{filename}")
    logger = MatcherLogger()
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filename)
    graph_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                logger=logger)
    start_alignment = time.time()
    graph_matcher.find_mapping()
    end_alignment = time.time()
    print(f"Matching time: {end_alignment - start_alignment} seconds")
    total_alignment_time += end_alignment - start_alignment
    count += 1
print(f"Average inference time {total_alignment_time / count} seconds")
# genealogical_graph = GenealogicalGraph(pedigree=pedigree)

# path_aware_genealogical_graph = PathAwareGenealogicalGraph(pedigree=pedigree)
# print(f"Number of levels: {len(path_aware_genealogical_graph.levels)}")
# descendant_cache = DescendantMemoryCache()
# descendant_cache = DescendantNoMemoryCache(path_aware_genealogical_graph)
# path_aware_genealogical_graph.set_descendant_writer(descendant_cache)

# path_aware_genealogical_graph.initialize_genealogical_graph_from_probands()
# path_processed_graph = PathProcessedGraph(path_aware_genealogical_graph)
# path_processed_graph.get_graph_from_genealogical_graph(path_aware_genealogical_graph)

# count = genealogical_graph.calculate_vertices_in_descendant_writer()
