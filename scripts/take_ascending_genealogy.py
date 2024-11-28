from graph.coalescent_tree import CoalescentTree
from graph.genealogical_graph import GenealogicalGraph
from scripts.utility import get_non_existing_path, read_integers_from_csv_file, get_filepath, \
    get_natural_number_input_in_bounds


def take_ascending_pedigree(pedigree_path: str, result_filepath: str, probands: [int]):
    pedigree = GenealogicalGraph.get_diploid_graph_from_file(filepath=pedigree_path)
    pedigree.save_ascending_genealogy_to_file(filepath=result_filepath, probands=probands)


def take_ascending_tree(tree_path: str, result_filepath: str, probands: [int]):
    # TODO: Test the whole implementation
    print("Parsing the tree")
    tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path, initialize_levels=True)
    # tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path, initialize_levels=False)
    # print("Reducing to ascending genealogy")
    # tree.reduce_to_ascending_genealogy(probands=probands, recalculate_levels=True)
    print("Removing the unary nodes")
    tree.remove_unary_nodes()
    tree.save_to_file(filepath=result_filepath)


def ask_and_parse_probands_from_file():
    probands_csv_filepath = get_filepath("Specify the path to the probands csv file:")
    return read_integers_from_csv_file(probands_csv_filepath)


def run_interactive_session():
    running_mode_prompt = ("Specify the running mode:\n"
                           "1) Take pedigree ascending genealogy\n"
                           "2) Take coalescent tree ascending genealogy\n")
    running_mode = get_natural_number_input_in_bounds(input_request=running_mode_prompt,
                                                      lower_bound=1,
                                                      upper_bound=2)
    input_filepath = get_filepath("Specify the path to the input file:")
    result_filepath = get_non_existing_path("Specify the result's filepath:")
    probands = ask_and_parse_probands_from_file()
    match running_mode:
        case 1:
            take_ascending_pedigree(pedigree_path=input_filepath, result_filepath=result_filepath, probands=probands)
        case 2:
            take_ascending_tree(tree_path=input_filepath, result_filepath=result_filepath, probands=probands)


if __name__ == '__main__':
    run_interactive_session()
