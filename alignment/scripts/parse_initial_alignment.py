from collections import defaultdict

from alignment.graph_matcher import validate_and_parse_yaml, InitialAssignment, PloidType
from scripts.utility import *


def validate_and_reduce_parsed_initial_assignments(parsed_assignments: dict[int, InitialAssignment]):
    # Without pedigree
    pedigree_dict = defaultdict(list)
    for assignment in parsed_assignments.values():
        pedigree_dict[assignment.pedigree_id].append(assignment.ploid_type)
    for pedigree_id, ploid_types in pedigree_dict.items():
        if len(ploid_types) > 2:
            raise ValueError(f"Too many ploid types for pedigree_id {pedigree_id}: {ploid_types}")

        for phased_ploid_type in [PloidType.Paternal, PloidType.Maternal]:
            if ploid_types.count(PloidType.Paternal) > 1:
                raise ValueError(f"Duplicate {phased_ploid_type} ploid type for pedigree_id {pedigree_id}")
        # TODO: Make pedigree_dict a matrix and use it to reduce the initial assignments


def main():
    initial_alignment_filepath = get_filepath("Specify the path to the initial assignments:")
    parsed_assignments = validate_and_parse_yaml(initial_alignment_filepath)
    try:
        validate_and_reduce_parsed_initial_assignments(parsed_assignments)
        print("Valid initial assignments")
        print(parsed_assignments)
    except ValueError as exception:
        print(f"Initial assignments are invalid: {exception}")


if __name__ == "__main__":
    main()
