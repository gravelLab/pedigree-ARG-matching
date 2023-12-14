import time

from genealogical_graph import *
from descendant_memory_cache import DescendantMemoryCache
from descendant_no_memory_cache import DescendantNoMemoryCache
from path_matching.path_processed_graph import PathAwareGenealogicalGraph

filename_input = input("Specify the pedigree file:")
filename = f"pedigrees/{filename_input}.pedigree"
pedigree_memory = Graph.get_pedigree_from_file(filename=filename)
pedigree_no_memory = Graph.get_pedigree_from_file(filename=filename)
genealogical_graph_memory = PathAwareGenealogicalGraph(pedigree=pedigree_memory)
genealogical_graph_no_memory = PathAwareGenealogicalGraph(pedigree=pedigree_no_memory)
descendant_cache_memory = DescendantMemoryCache()
descendant_cache_no_memory = DescendantNoMemoryCache(genealogical_graph_no_memory)
genealogical_graph_memory.set_descendant_writer(descendant_cache_memory)
genealogical_graph_no_memory.set_descendant_writer(descendant_cache_no_memory)
start = time.time()
genealogical_graph_memory.initialize_genealogical_graph_from_probands()
end = time.time()
print(f"Running time (memory): {end - start} seconds")

start = time.time()
genealogical_graph_no_memory.initialize_genealogical_graph_from_probands()
end = time.time()
print(f"Running time (no memory): {end - start} seconds")

