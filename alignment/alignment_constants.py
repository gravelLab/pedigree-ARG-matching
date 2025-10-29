from enum import Enum

from alignment.configuration import AlignmentVertexMode, AlignmentEdgeMode, PosteriorProbabilitiesCalculationMode


class PloidType(Enum):
    Paternal = "P"
    Maternal = "M"


ploid_types = (PloidType.Paternal.value, PloidType.Maternal.value)


# YAML keys used in the driver file
initial_assignments_key = "initial_assignments"
coalescent_tree_key = "coalescent_tree"
pedigree_key = "pedigree"
coalescent_id_key = "coalescent_id"
path_key = "path"
pedigree_ids_key = "pedigree_ids"
separation_symbol_key = "separation_symbol"
missing_parent_notation_key = "missing_parent_notation"
skip_first_line_key = "skip_first_line"
verify_graph_has_no_cycles_key = "check_for_cycles"
output_path_key = "output_path"
alignment_vertex_mode_key = "alignment_vertex_mode"
alignment_edge_mode_key = "alignment_edge_mode"
posterior_probability_calculation_mode_key = "posterior_probability_calculation_mode"

alignment_vertex_mode_dict = {
    "all": AlignmentVertexMode.ALL_ALIGNMENTS,
    "example_per_root_assignment": AlignmentVertexMode.EXAMPLE_PER_ROOT_ASSIGNMENT
}

alignment_edge_mode_dict = {
    "one": AlignmentEdgeMode.EXAMPLE_EDGE_ALIGNMENT,
    "all": AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS
}

posterior_probability_calculation_mode_dict = {
    "skip": PosteriorProbabilitiesCalculationMode.SKIP,
    "vertex_alignment": PosteriorProbabilitiesCalculationMode.VERTEX_ALIGNMENT_PROBABILITY,
    "vertex_inclusion": PosteriorProbabilitiesCalculationMode.INDIVIDUAL_INCLUSION_PROBABILITY
}
