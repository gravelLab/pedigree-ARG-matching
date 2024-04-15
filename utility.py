import os
import matplotlib.pyplot as plt


def get_file_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif not os.path.isfile(file_path):
            print("The specified path is not a file, try again")
        else:
            return file_path


def get_directory_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif os.path.isfile(file_path):
            print("The specified path is a file, try again")
        else:
            return file_path


def get_non_existing_directory_name(input_request: str):
    while True:
        directory_name = input(input_request)
        if os.path.exists(directory_name):
            print("The specified directory already exists, try again")
        else:
            return directory_name


def get_integer_input(input_request: str):
    while True:
        try:
            result = int(input(input_request))
            if result < 1:
                print("Specify a positive value")
            else:
                return result
        except ValueError:
            print("You need to specify an integer")


def parse_dictionary_from_file(file_path: str):
    result = dict()
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            key = int(parts[0].strip())
            value = int(parts[1].strip())
            result[key] = value
    return result


def build_histogram(histogram_filename: str, dictionary: dict):
    x = list(dictionary.keys())
    y = list(dictionary.values())

    plt.bar(x, y)
    plt.xlabel('Distances')
    plt.ylabel('Frequency')
    plt.title('Distance histogram')
    plt.savefig(f"{histogram_filename}.png")
    plt.close()
