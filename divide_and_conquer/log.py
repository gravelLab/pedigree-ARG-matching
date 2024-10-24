from divide_and_conquer.configuration import print_enabled


def print_log(log_string: str):
    if print_enabled:
        print(log_string)
