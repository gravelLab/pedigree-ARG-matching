from enum import Enum

section_separator = 60 * "#" + '\n'
subsection_separator = 60 * "-" + '\n'


class ProbandInitialAssignmentsMode(Enum):
    PLOID = 1,
    INDIVIDUAL = 2,


class AlignmentVertexMode(Enum):
    ALL_ALIGNMENTS = 1,
    EXAMPLE_PER_ROOT_ASSIGNMENT = 2,


class AlignmentEdgeMode(Enum):
    EXAMPLE_EDGE_ALIGNMENT = 1,
    ALL_EDGE_ALIGNMENTS = 2,


default_proband_initial_assignments_mode = ProbandInitialAssignmentsMode.INDIVIDUAL
default_alignment_vertex_mode = AlignmentVertexMode.ALL_ALIGNMENTS
default_alignment_edge_mode = AlignmentEdgeMode.ALL_EDGE_ALIGNMENTS

# Alignment statistics section start
calculate_distances_histogram = False
default_calculate_posterior_probabilities = True
calculate_similarity = False
save_edge_alignments = False
# Alignment statistics section end

default_missing_parent_notation = "-1"
default_separation_symbol = " "
default_skip_first_line = False

logs_enabled = True
logs_default_directory_name = "logs"

# File extensions
pedigree_extension = ".pedigree"
tree_extension = ".tree"

print_enabled = True
