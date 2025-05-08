from alignment.driver_file import ParsedDriverFile
from scripts.utility import *


def main():
    initial_alignment_filepath = get_filepath("Specify the path to the initial assignments:")
    parsed_assignments = ParsedDriverFile.parse_driver_file_and_validate_initial_assignments(initial_alignment_filepath)
    try:
        print("Valid initial assignments")
        print(parsed_assignments)
    except ValueError as exception:
        print(f"The driver file is invalid: {exception}")


if __name__ == "__main__":
    main()
