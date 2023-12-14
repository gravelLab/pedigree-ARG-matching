from path_matching.path_processed_graph import PathAwareGenealogicalGraph


class Variable:
    def __init__(self, proband_pedigree: int, ancestor_coalescent_tree: int,
                 domain):
        self.proband_pedigree = proband_pedigree
        self.ancestor_coalescent_tree = ancestor_coalescent_tree
        self.domain = domain
        self.value = None

    def get_pedigree_candidates(self):
        return self.domain.keys()
