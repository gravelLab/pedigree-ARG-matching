from alignment.graph_matcher import validate_and_parse_yaml
from scripts.utility import *


def main():
    initial_alignment_filepath = get_filepath("Specify the path to the initial assignments:")
    parsed_assignments = validate_and_parse_yaml(initial_alignment_filepath)
    try:
        print("Valid initial assignments")
        print(parsed_assignments)
    except ValueError as exception:
        print(f"Initial assignments are invalid: {exception}")


if __name__ == "__main__":
    main()
