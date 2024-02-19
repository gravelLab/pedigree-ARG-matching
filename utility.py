import os


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
