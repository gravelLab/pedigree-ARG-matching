from enum import Enum


class InitialMatchingMode(Enum):
    PLOID = 1
    INDIVIDUAL = 2


class MatchingMode(Enum):
    ALL_ALIGNMENTS = 1,
    EXAMPLE_PER_ROOT_ASSIGNMENT = 2,


default_initial_matching_mode = InitialMatchingMode.INDIVIDUAL
current_matching_mode = MatchingMode.EXAMPLE_PER_ROOT_ASSIGNMENT

logs_enabled = True
logs_default_directory_name = "logs"

print_enabled = True
