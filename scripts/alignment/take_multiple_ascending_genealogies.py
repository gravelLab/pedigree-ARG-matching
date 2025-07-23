import warnings
from concurrent.futures import ProcessPoolExecutor

from lineagekit.core.CoalescentTree import CoalescentTree
from lineagekit.core.PloidPedigree import PloidPedigree

from scripts.utility.basic_utility import *


def take_ascending_genealogies(input_genealogies_directory_path: str, probands: [int], result_directory_filepath: str,
                               processing_function, processes_number: int = 1):
    input_genealogies_directory_path = Path(input_genealogies_directory_path)
    result_directory_filepath = Path(result_directory_filepath)
    result_directory_filepath.mkdir(parents=True)
    input_genealogy_entries = os.listdir(input_genealogies_directory_path)
    with ProcessPoolExecutor(max_workers=processes_number) as executor:
        futures = [
            executor.submit(
                processing_function,
                input_genealogies_directory_path,
                result_directory_filepath,
                genealogy_entry,
                probands,
            )
            for genealogy_entry in input_genealogy_entries
        ]
        for future in futures:
            future.result()


def process_pedigree_directory(parent_pedigrees_directory_path: Path, result_directory_filepath: Path,
                               pedigree_directory: str, probands: [int]):
    # Paths for the current pedigree directory
    pedigree_directory_path = parent_pedigrees_directory_path / pedigree_directory
    result_pedigree_directory_path = result_directory_filepath / pedigree_directory
    os.makedirs(result_pedigree_directory_path, exist_ok=True)
    # Find and process the pedigree file
    try:
        pedigree_filename = get_unique_filename_with_specified_extension(pedigree_directory_path, pedigree_extension)
    except ValueError as ex:
        warnings.warn(ex)
        return
    pedigree_filepath = pedigree_directory_path / pedigree_filename
    print(f"Processing {pedigree_filepath}")
    pedigree = PloidPedigree.get_ploid_pedigree_from_file(filepath=pedigree_filepath)
    ascending_pedigree_filepath = result_pedigree_directory_path / pedigree_filename
    pedigree.save_ascending_genealogy_as_diploid(filepath=ascending_pedigree_filepath, vertices=probands)


def process_coalescent_tree_directory(parent_coalescent_tree_directory_path: Path, result_directory_filepath: Path,
                                      tree_filename: str, probands: [int]):
    input_tree_filepath = parent_coalescent_tree_directory_path / tree_filename
    result_tree_filepath = result_directory_filepath / tree_filename
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=input_tree_filepath)
    os.makedirs(result_directory_filepath, exist_ok=True)
    coalescent_tree.save_ascending_graph_to_file(filepath=result_tree_filepath, vertices=probands)


def process_pedigree_tree_directory(pedigree_tree_parent_directory_path: Path, result_directory_filepath: Path,
                                    tree_pedigree_directory: str, probands: [int]):
    tree_pedigree_directory_path = pedigree_tree_parent_directory_path / tree_pedigree_directory
    process_pedigree_directory(pedigree_directory=tree_pedigree_directory,
                               result_directory_filepath=result_directory_filepath,
                               probands=probands,
                               parent_pedigrees_directory_path=pedigree_tree_parent_directory_path)
    non_pedigree_files = [f for f in os.listdir(tree_pedigree_directory_path) if not f.endswith(pedigree_extension)]
    # Handle warnings
    if len(non_pedigree_files) == 0:
        warnings.warn(f"No non-pedigree files found in the {tree_pedigree_directory_path} directory")
        return
    if len(non_pedigree_files) > 1:
        warnings.warn(f"Multiple non-pedigree files found in the {tree_pedigree_directory_path} directory")
        return
    coalescent_tree_filename = non_pedigree_files[0]
    # Process the input files
    process_coalescent_tree_directory(parent_coalescent_tree_directory_path=tree_pedigree_directory_path,
                                      result_directory_filepath=result_directory_filepath,
                                      tree_filename=coalescent_tree_filename,
                                      probands=probands)


def ask_and_parse_probands_from_file():
    probands_csv_filepath = get_filepath("Specify the path to the probands csv file:")
    return read_integers_from_csv_file(probands_csv_filepath)


def run_interactive_session_pedigrees_ascending():
    pedigrees_directory_path = get_directory_path("Specify the path to the pedigree directory:")
    result_directory_filepath = get_non_existing_path("Specify the path where the results will be saved:")
    probands = ask_and_parse_probands_from_file()
    processes_number = get_natural_number_input("Specify the number of processes to use:")
    take_ascending_genealogies(input_genealogies_directory_path=pedigrees_directory_path,
                               probands=probands,
                               result_directory_filepath=result_directory_filepath,
                               processing_function=process_pedigree_directory,
                               processes_number=processes_number)


def run_interactive_session_trees_ascending():
    trees_directory_path = get_directory_path("Specify the path to the tree directory:")
    result_directory_filepath = get_non_existing_path("Specify the path where the results will be saved:")
    probands = ask_and_parse_probands_from_file()
    processes_number = get_natural_number_input("Specify the number of processes to use:")
    take_ascending_genealogies(input_genealogies_directory_path=trees_directory_path,
                               result_directory_filepath=result_directory_filepath,
                               processing_function=process_coalescent_tree_directory,
                               probands=probands, processes_number=processes_number)


def run_interactive_session_pedigree_tree_directories_ascending():
    tree_pedigree_parent_directory_path = get_directory_path("Specify the path to the parent directory of "
                                                             "pedigree-tree directories:")
    result_directory_filepath = get_non_existing_path("Specify the path where the results will be saved:")
    probands = ask_and_parse_probands_from_file()
    processes_number = get_natural_number_input("Specify the number of processes to use:")
    take_ascending_genealogies(input_genealogies_directory_path=tree_pedigree_parent_directory_path,
                               result_directory_filepath=result_directory_filepath,
                               processing_function=process_pedigree_tree_directory,
                               probands=probands, processes_number=processes_number)


def run_interactive_session():
    input_prompt = ("Specify the running mode:\n"
                    "1) Take ascending genealogies for pedigree directories\n"
                    "2) Take ascending genealogies for coalescent trees\n"
                    "3) Take ascending genealogies for pedigree-tree directories\n")
    selected_mode = get_natural_number_input_in_bounds(input_prompt, 1, 3)
    match selected_mode:
        case 1:
            run_interactive_session_pedigrees_ascending()
        case 2:
            run_interactive_session_trees_ascending()
        case 3:
            run_interactive_session_pedigree_tree_directories_ascending()


if __name__ == '__main__':
    run_interactive_session()
