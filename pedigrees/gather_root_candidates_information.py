import os

from utility import build_histogram, parse_dictionary_from_file

root_directory = "11_19_simulation"
resulting_statistics_filename = f"{root_directory}_root_candidates_number.txt"
test_subdirectory = "test"


def get_number_of_candidates(file_path: str, clade_root: str):
    if not os.path.exists(file_path):
        return None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(f"{clade_root}:"):
                # Parsing the remaining text after "49504682:"
                parsed_text = int(line.split(f"{clade_root}:", 1)[1].strip())
                return parsed_text
    return None


def gather_root_candidates_number(root_dir: str, result_filename: str):
    with open(f"{result_filename}", "w") as file:
        for directory in os.listdir(root_dir):
            if not os.path.isdir(os.path.join(root_dir, directory)):
                continue
            for subdirectory in os.listdir(os.path.join(root_dir, directory)):
                clade_directory_path = os.path.join(root_dir, directory, subdirectory, test_subdirectory, "clade")
                for clade_directory in os.listdir(clade_directory_path):
                    if not os.path.isdir(os.path.join(clade_directory_path, clade_directory)):
                        continue
                    if clade_directory == "logs":
                        continue
                    target_file = f"_clade_{clade_directory}.txt"
                    result_filepath = os.path.join(clade_directory_path,
                                                   clade_directory, target_file)
                    number_of_candidates = get_number_of_candidates(result_filepath, clade_directory)
                    if number_of_candidates is not None:
                        file.write(f"{directory}: {number_of_candidates}\n")


def build_histograms(root_dir: str):
    for directory in os.listdir(root_dir):
        if not os.path.isdir(os.path.join(root_dir, directory)):
            continue
        for subdirectory in os.listdir(os.path.join(root_dir, directory)):
            clade_directory_path = os.path.join(root_dir, directory, subdirectory, test_subdirectory, "clade")
            if not os.path.exists(clade_directory_path):
                continue
            for clade_directory in os.listdir(clade_directory_path):
                if not os.path.isdir(os.path.join(clade_directory_path, clade_directory)):
                    continue
                if clade_directory == "logs":
                    continue
                histogram_filename = f"distance_histogram_{clade_directory}.txt"
                histogram_filepath = os.path.join(clade_directory_path, clade_directory, histogram_filename)
                if not os.path.exists(histogram_filepath):
                    continue
                histogram_dictionary = parse_dictionary_from_file(histogram_filepath)
                histogram_image_filepath = os.path.join(clade_directory_path, clade_directory,
                                                        f"distance_histogram_{clade_directory}")
                build_histogram(histogram_filename=histogram_image_filepath, dictionary=histogram_dictionary)


# build_histograms(root_directory)
gather_root_candidates_number(root_directory, resulting_statistics_filename)
