from abc import abstractmethod
from genealogical_graph import CoalescentTree, GenealogicalGraph


class CommonAncestorQueryableGraph:

    @abstractmethod
    def get_common_ancestors(self, first_vertex: int, second_vertex: int):
        pass

    @abstractmethod
    def verify_alignment(self, coalescent_tree: CoalescentTree, pedigree: GenealogicalGraph, mapping: dict):
        pass
