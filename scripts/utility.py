import os
from pathlib import Path

import matplotlib.pyplot as plt
from scipy.stats import poisson
import random
from math import ceil


def get_filepaths(input_request: str) -> list[str]:
    print(input_request)
    paths = []
    while True:
        path = input("Enter a file path (or press Enter to finish): ").strip()
        if path == "":
            if paths:
                break
            else:
                print("Please enter at least one file path.")
        elif os.path.isfile(path):
            paths.append(path)
        else:
            print("Invalid file path. Please try again.")
    return paths


def get_directory_paths(input_request: str) -> list[str]:
    print(input_request)
    paths = []
    while True:
        path = input("Enter a folder path (or press Enter to finish): ").strip()
        if path == "":
            if paths:
                break
            else:
                print("Please enter at least one folder path.")
        elif os.path.isdir(path):
            paths.append(path)
        else:
            print("Invalid folder path. Please try again.")
    return paths


def get_filepath(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif not os.path.isfile(file_path):
            print("The specified path is not a file, try again")
        else:
            return file_path


def verify_filepath(path: str) -> bool:
    return os.path.exists(path) and os.path.isfile(path)


def verify_folder_path(path: str) -> bool:
    return os.path.exists(path) and os.path.isdir(path)


def get_directory_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif os.path.isfile(file_path):
            print("The specified path is a file, try again")
        else:
            return file_path


def get_non_empty_string(input_request: str) -> str:
    while True:
        response = input(input_request)
        if not response:
            print("Specify a non-empty string")
        else:
            return response


def get_non_existing_path(input_request: str):
    while True:
        path = input(input_request)
        if not path:
            print("Specify a non-empty string")
        elif os.path.exists(path):
            print("The specified path exists, try again")
        else:
            return path


def get_natural_number_input(input_request: str):
    while True:
        try:
            result = int(input(input_request))
            if result < 1:
                print("Specify a positive value")
            else:
                return result
        except ValueError:
            print("You need to specify an integer")


def get_natural_number_input_in_bounds(input_request: str, lower_bound: int, upper_bound: int):
    if lower_bound >= upper_bound:
        raise ValueError("Lower bound cannot be greater than upper bound")
    while True:
        value = get_natural_number_input(input_request)
        if value < lower_bound or value > upper_bound:
            print("Value out of bounds, try again")
            continue
        return value


def get_natural_number_with_lower_bound(input_request: str, lower_bound: int):
    while True:
        value = get_natural_number_input(input_request)
        if value < lower_bound:
            print("Value out of bounds, try again")
            continue
        return value


def save_dictionary_to_file(dictionary_filepath: str, dictionary: dict):
    dictionary_file = open(dictionary_filepath, 'w')
    for key, value in dictionary.items():
        dictionary_file.write(f"{key}: {value}\n")
    dictionary_file.close()


def parse_dictionary_from_file(file_path: str):
    result = dict()
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            key = int(parts[0].strip())
            value = int(parts[1].strip())
            result[key] = value
    return result


def read_mapping_from_file(filepath):
    parsed_dict = {}
    with open(filepath, 'r') as file:
        for line in file:
            # Split the line on ':' to separate key and value
            key, value = line.strip().split(':')
            # Convert key to an integer and value to a list of integers
            key = int(key.strip())
            value = [int(x) for x in value.strip()[1:-1].split(',')]
            # Add to the dictionary
            parsed_dict[key] = value
    return parsed_dict


def build_histogram(histogram_filepath: str, dictionary: dict):
    x = list(dictionary.keys())
    y = list(dictionary.values())

    plt.bar(x, y)
    plt.xlabel('Distances')
    plt.ylabel('Frequency')
    plt.title('Distance histogram')
    plt.savefig(f"{histogram_filepath}.png")
    plt.close()


def random_subselect(input_list, percentage):
    # Calculate the number of elements to select
    num_elements = ceil(len(input_list) * percentage)
    # Randomly select elements from the list
    return random.sample(input_list, num_elements)


def random_subselect_poisson(input_list, percentage):
    n = len(input_list)
    math_expectation = n * percentage
    number_of_errors = poisson.rvs(math_expectation)
    return random.sample(input_list, number_of_errors)


def get_yes_or_no(prompt: str) -> bool:
    while True:
        user_input = input(f"{prompt} (yes or no): ").strip().lower()
        if user_input in ['yes', 'y']:
            return True
        elif user_input in ['no', 'n']:
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")


def read_integers_from_csv_file(filepath: str):
    with open(filepath, mode="r") as file:
        content = file.read().strip()
        integer_list = [int(value) for value in content.split(",")]
    return integer_list


def get_filenames_with_given_extension(directory_path: str, extension: str):
    return [file for file in os.listdir(directory_path) if file.endswith(extension)]


def get_unique_filename_with_specified_extension(directory_path: str, extension: str):
    files = get_filenames_with_given_extension(directory_path=directory_path, extension=extension)
    if len(files) == 0:
        raise Exception(f"There are no {extension} files in the directory")
    if len(files) == 2:
        raise Exception(f"There are multiple {extension} files in the directory")
    return files[0]


def dict_has_duplicate_values(dictionary: dict):
    seen = set()
    for value in dictionary.values():
        if value in seen:
            return True
        seen.add(value)
    return False


def get_paths_from_tree_pedigree_directory(tree_pedigree_directory_path: str | Path):
    tree_pedigree_directory_path = Path(tree_pedigree_directory_path)
    if not tree_pedigree_directory_path.is_dir():
        return None
    files = list(tree_pedigree_directory_path.glob("*"))
    files = [file.resolve() for file in files if file.is_file()]
    if len(files) != 2:
        return None
    pedigree_filename = next((file for file in files if file.suffix == ".pedigree"), None)
    tree_filename = next((file for file in files if file != pedigree_filename), None)
    pedigree_filepath = tree_pedigree_directory_path / pedigree_filename
    tree_filepath = tree_pedigree_directory_path / tree_filename
    return pedigree_filepath, tree_filepath
