from io import UnsupportedOperation
import traceback
from divide_and_conquer.graph_matcher import *
from run_alignment import save_alignment_result_to_files, get_alignments_proband_distance_probands_ignored
from utility import get_file_path, get_directory_path, get_non_existing_directory_name, os

log_directory = "log_dir"
initial_alignment_dir_name = "initial_alignment"
general_result_filename = "!simulation_result.txt"
detailed_result_filename = "!detailed_simulation_result.txt"

# Keeping track of the number of simulations with 0 solutions
no_solutions_counter = 0
# Keeping track of the number of simulations where the root assignments are the root vertex itself and their spouses
# In other words, it is a subset of the correct assignments
only_individual_and_spouses_counter = 0
# The number of simulations where the correct individual is found, but a non-spouse assignment
# for the root is also present
individual_and_non_spouse_present_counter = 0
# The number of simulations where the individual is not present, but one of the spouses is present
individual_not_present_spouse_present_counter = 0
# The number of simulations where neither the individual nor a spouse is present
neither_individual_nor_spouse_present_counter = 0

# Other general simulation counters
total_number_of_solutions = 0
total_distance_to_solution = 0
total_number_of_root_vertex_individual_assignments = 0
new_individual_assignments_number = 0
new_individual_simulations_number = 0


def log_graph_paths(pedigree_path: str, tree_path: str):
    general_result_file.write(f"The path to the initial, error-free pedigree is: {pedigree_path}\n")
    general_result_file.write(f"The path to the coalescent tree is: {tree_path}\n")
    general_result_file.write("#####################################\n")


def log_solutions(alignments: list[dict], simulation_step: int):
    detailed_result_file.write("#####################################\n")
    detailed_result_file.write(f"Simulation {simulation_step}\n")
    detailed_result_file.write(f"Number of incorrect alignments in this simulation: {len(alignments)}\n")
    for incorrect_alignment in alignments:
        detailed_result_file.write("------------------------------------\n")
        detailed_result_file.write(f"{incorrect_alignment}\n")
        detailed_result_file.write("------------------------------------\n")
    detailed_result_file.write("#####################################\n")
    detailed_result_file.flush()


def log_overall_results(number_of_simulation_steps: int):
    if total_number_of_solutions == 0:
        # Nothing to log
        return
    simulations_with_solutions = number_of_simulation_steps - no_solutions_counter
    alignments_to_identity_average_distance = (total_distance_to_solution / total_number_of_solutions)
    delta = alignments_to_identity_average_distance - initial_average_distance_to_identity
    number_of_probands = len(coalescent_tree.probands)
    non_proband_vertices = len(coalescent_tree.vertex_to_level_map) - number_of_probands
    general_result_file.write(
        f"Initial average distance to the identity is: {initial_average_distance_to_identity}\n")
    general_result_file.write(f"Average distance from simulation solutions to the identity "
                              f"is: {alignments_to_identity_average_distance}\n")
    general_result_file.write(f"Average distance delta is: {delta}\n")
    general_result_file.write(f"The number of probands is {number_of_probands}\n")
    general_result_file.write(f"The number of non-proband vertices is: {non_proband_vertices}\n")
    general_result_file.write(
        f"Average distance is {alignments_to_identity_average_distance / non_proband_vertices} "
        f"from the number of non-proband vertices in the coalescent tree\n")
    general_result_file.write("----------------------------------------------------------------------\n")
    general_result_file.write(f"The number of simulations with solutions: {simulations_with_solutions}\n")
    general_result_file.write(f"The number of simulations with no solutions: {no_solutions_counter}\n")
    general_result_file.write(f"The number of simulations with only the individual and spouses: "
                              f"{only_individual_and_spouses_counter}\n")
    general_result_file.write(f"The number of simulations with the individual and a non-spouse present: "
                              f"{individual_and_non_spouse_present_counter}\n")
    general_result_file.write(f"The number of simulations without the individual, but with a spouse: "
                              f"{individual_not_present_spouse_present_counter}\n")
    general_result_file.write(f"The number of simulations without the individual and all the spouses:"
                              f"{neither_individual_nor_spouse_present_counter} ")
    general_result_file.write("----------------------------------------------------------------------\n")
    general_result_file.write(f"Initial number of individual candidates for the root: "
                              f"{initial_root_vertex_individual_assignments_number}\n")
    general_result_file.write(
        f"Average number of individual candidates for the root per simulation with solutions: "
        f"{total_number_of_root_vertex_individual_assignments / simulations_with_solutions}\n")
    general_result_file.write(f"New individual root assignments: {new_individual_assignments_number}\n")
    general_result_file.write(f"Average number of new individual root assignments per simulation with solutions: "
                              f"{new_individual_assignments_number / simulations_with_solutions}\n")
    general_result_file.write(f"Number of new individual root assignments simulations: "
                              f"{new_individual_simulations_number}\n")


coalescent_tree_path = get_file_path("Specify the absolute path to the coalescent tree:")
coalescent_tree_filename = os.path.basename(coalescent_tree_path)
pedigree_no_errors_path = get_file_path("Specify the absolute path to the initial, error-free pedigree:")

error_pedigrees_folder = get_directory_path("Specify the absolute path to the error pedigrees folder:")
os.chdir(error_pedigrees_folder)
os.chdir("..")
simulation_name = get_non_existing_directory_name("Specify the simulation name:")
os.mkdir(simulation_name)
os.chdir(simulation_name)
detailed_result_file = open(detailed_result_filename, 'w')
general_result_file = open(general_result_filename, 'w')
log_graph_paths(pedigree_path=pedigree_no_errors_path, tree_path=coalescent_tree_path)

coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filename=coalescent_tree_path)
coalescent_tree.remove_unary_nodes()

# Run the initial alignment
initial_pedigree = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filename=pedigree_no_errors_path,
                                                                             initialize_levels=False,
                                                                             initialize_ancestor_maps=False
                                                                             )
if current_initial_matching_mode == InitialMatchingMode.PLOID:
    ascending_genealogy_probands = coalescent_tree.probands
else:
    proband_individuals = {proband // 2 for proband in coalescent_tree.probands}
    probands_all_ploids = [
        ploid
        for proband_individual in proband_individuals
        for ploid in [2 * proband_individual, 2 * proband_individual + 1]
    ]
    ascending_genealogy_probands = probands_all_ploids
initial_pedigree.reduce_to_ascending_genealogy(probands=ascending_genealogy_probands,
                                               recalculate_levels=True)
initial_pedigree.initialize_potential_mrca_map()
os.mkdir(initial_alignment_dir_name)
os.chdir(initial_alignment_dir_name)
initial_logger = MatcherLogger(logs_directory_name=log_directory)
initial_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=initial_pedigree,
                              logger=initial_logger)
initial_result_dict = initial_matcher.find_mapping()
save_alignment_result_to_files(initial_result_dict, coalescent_tree, initial_pedigree, coalescent_tree_filename)

os.chdir("..")
number_of_clades = len(initial_result_dict)
if number_of_clades != 1:
    raise UnsupportedOperation("The current version of the script supports trees with only one clade")
print(f"Clades: {number_of_clades}")
root_vertex = next(iter(initial_result_dict))
root_vertex_individual = root_vertex // 2
root_vertex_children_all_ploids = [
    child_ploid
    for child in initial_pedigree.children_map[root_vertex]
    for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
]
# In the phased case, we cannot guarantee that the other ploid of every root vertex child ploid is present in the graph
individual_and_spouses = {p // 2 for child_ploid in root_vertex_children_all_ploids
                          for p in initial_pedigree.parents_map.get(child_ploid, [])}
assert individual_and_spouses, "Zero assignments found for the root!"
valid_alignments = [d for lst in initial_result_dict.values() for d in lst]
initial_result_dict = None
initial_pedigree = None
# Gathering the statistics about the original alignments
initial_root_vertex_individual_assignments_number = 0
identity_solution = coalescent_tree.get_identity_solution()
assert identity_solution in valid_alignments
initial_average_distance_to_identity = sum(
    get_alignments_proband_distance_probands_ignored(alignment, identity_solution, coalescent_tree.probands)
    for alignment in valid_alignments
) / len(valid_alignments)
general_result_file.write(f"Total number of initial alignments is: {len(valid_alignments)}\n")
general_result_file.write(f"Average distance to the identity is: {initial_average_distance_to_identity}\n")
# Calculating the number of individuals to which the clade's root is assigned
initial_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in valid_alignments}
initial_root_vertex_individual_assignments_number += len(initial_root_vertex_individual_assignments)
general_result_file.write(f"Number of individual assignments for the root:"
                          f" {initial_root_vertex_individual_assignments_number}\n")
general_result_file.write(f"Clade root individual: {root_vertex_individual}\n")
general_result_file.write(f"Clade root and their spouses: {individual_and_spouses}\n")
general_result_file.write(f"Individual assignments for the root: {initial_root_vertex_individual_assignments}\n")
general_result_file.flush()

# Run the alignments on the pedigrees with errors
os.chdir(error_pedigrees_folder)
current_directory = os.getcwd()

print("Starting running the alignments")
simulation_counter = 0
try:
    for subdirectory in os.listdir(current_directory):
        os.chdir(subdirectory)
        simulation_counter += 1
        for file in os.listdir(os.getcwd()):
            if file.endswith(".pedigree"):
                print(f"Parsing {file}")
                potential_mrca_graph = (PotentialMrcaProcessedGraph.
                                        get_processed_graph_from_file(filename=file,
                                                                      initialize_levels=False,
                                                                      initialize_ancestor_maps=False))
                print("Reducing to the ascending genealogy")
                potential_mrca_graph.reduce_to_ascending_genealogy(probands=ascending_genealogy_probands,
                                                                   recalculate_levels=True)
                potential_mrca_graph.initialize_potential_mrca_map()
                os.chdir(f"../../{simulation_name}")
                os.mkdir(subdirectory)
                os.chdir(subdirectory)
                logger = MatcherLogger(logs_directory_name=log_directory)
                matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                      logger=logger)
                result = matcher.find_mapping()
                save_alignment_result_to_files(result, coalescent_tree, potential_mrca_graph, coalescent_tree_filename)
                resulting_alignments = [d for lst in result.values() for d in lst]
                log_solutions(resulting_alignments, simulation_counter)
                general_result_file.write("#####################################\n")
                general_result_file.write(f"Simulation step: {simulation_counter}\n")
                general_result_file.write(f"Number of alignments in this simulation {len(resulting_alignments)}\n")
                if len(resulting_alignments) > 0:
                    total_number_of_solutions += len(resulting_alignments)
                    total_distance_to_solution += sum(
                        get_alignments_proband_distance_probands_ignored(alignment, identity_solution,
                                                                         coalescent_tree.probands)
                        for alignment in resulting_alignments
                    )
                    general_result_file.write(f"New average distance from simulation solutions to the"
                                              f" identity is: "
                                              f"{total_distance_to_solution / total_number_of_solutions}\n")
                    # Calculate the number of root assignments for this simulation
                    simulation_root_vertex_individual_assignments = {x[root_vertex] // 2 for x in resulting_alignments}
                    general_result_file.write(f"Root vertex individual assignments "
                                              f"({len(simulation_root_vertex_individual_assignments)}): "
                                              f"{simulation_root_vertex_individual_assignments}\n")
                    total_number_of_root_vertex_individual_assignments += (
                        len(simulation_root_vertex_individual_assignments))
                    new_individual_assignments_number += len({x for x in simulation_root_vertex_individual_assignments
                                                              if x not in initial_root_vertex_individual_assignments})
                    if new_individual_assignments_number > 0:
                        new_individual_simulations_number += 1
                    # Classify the simulation result
                    if root_vertex_individual in simulation_root_vertex_individual_assignments:
                        if simulation_root_vertex_individual_assignments.issubset(individual_and_spouses):
                            only_individual_and_spouses_counter += 1
                        else:
                            individual_and_non_spouse_present_counter += 1
                    else:
                        if simulation_root_vertex_individual_assignments.intersection(individual_and_spouses):
                            individual_not_present_spouse_present_counter += 1
                        else:
                            neither_individual_nor_spouse_present_counter += 1
                else:
                    no_solutions_counter += 1
                general_result_file.write("#####################################\n")
                general_result_file.flush()
                os.chdir(current_directory)
                break
except KeyboardInterrupt:
    print("Stop the simulation, log the final results")
except Exception as ex:
    print("Exception occurred:")
    traceback.print_exc()

log_overall_results(number_of_simulation_steps=simulation_counter)
