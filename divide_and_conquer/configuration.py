from enum import Enum


class InitialMatchingMode(Enum):
    PLOID = 1
    INDIVIDUAL = 2


class MatchingMode(Enum):
    ALL_ALIGNMENTS = 1,
    LIKELIHOOD = 2,


current_initial_matching_mode = InitialMatchingMode.INDIVIDUAL
current_matching_mode = MatchingMode.ALL_ALIGNMENTS

logs_enabled = True
logs_default_directory_name = "logs"

print_enabled = True
