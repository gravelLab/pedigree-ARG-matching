from enum import Enum


class InitialMatchingMode(Enum):
    PLOID = 1
    INDIVIDUAL = 2


class MatchingMode(Enum):
    ALL_ALIGNMENTS = 1,
    EXAMPLE_PER_ROOT_ASSIGNMENT = 2,


default_initial_mapping_mode = InitialMatchingMode.INDIVIDUAL
default_alignment_mode = MatchingMode.EXAMPLE_PER_ROOT_ASSIGNMENT

# Alignment statistics section start
calculate_distances_histogram = False
calculate_likelihood = False
calculate_similarity = False
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
