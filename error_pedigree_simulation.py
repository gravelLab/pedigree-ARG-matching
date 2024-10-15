import copy
import random
import warnings
from math import ceil

from divide_and_conquer.graph_matcher import *
from run_alignment import save_alignment_result_to_files, get_alignments_proband_distance_probands_ignored
from utility import *

from scipy.stats import poisson

simulation_root_dir = "error_pedigree_simulation"
log_directory = "log_dir"


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


error_rate = 0.01

os.makedirs(simulation_root_dir, exist_ok=True)
os.chdir(simulation_root_dir)

simulation_name = get_non_existing_directory_name("Specify the name of the simulation: ")
simulation_steps = get_integer_input("Specify the number of simulation steps: ")
pedigree_path = get_file_path("Specify the absolute path to the pedigree:")
tree_path = get_file_path("Specify the absolute path to the coalescent tree:")
tree_filename = os.path.basename(tree_path)

os.makedirs(simulation_name)
os.chdir(simulation_name)

detailed_result_filename = "_detailed_simulation_result.txt"
detailed_result_file = open(detailed_result_filename, 'w')

# total_number_of_incorrect_solutions = 0
total_number_of_solutions = 0
total_distance_to_solution = 0
general_result_filename = "_simulation_result.txt"
general_result_file = open(general_result_filename, 'w')
# Counters
# Keeping track of the number of simulations with 0 solutions
no_solutions_counter = 0
# Keeping track of the number of simulations where the root assignments are the root vertex itself and their spouses
# In other words, it is a subset of the correct assignments
only_correct_root_assignments_counter = 0
# Keeping track of the number of simulations where the root assignments contain the root vertex itself or one of
# their spouses
present_correct_root_assignments_counter = 0
# Keeping track of the number of simulations where the root assignments don't contain the root vertex and
# all of their spouses
no_correct_root_assignments_counter = 0

general_result_file.write(f"The path to the pedigree is: {pedigree_path}\n")
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


potential_mrca_graph = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filename=pedigree_path,
                                                                                 initialize_ancestor_maps=False)
initial_parents_map = copy.deepcopy(potential_mrca_graph.parents_map)
initial_children_map = copy.deepcopy(potential_mrca_graph.children_map)

# min_levels, min_levels_dict = potential_mrca_graph.get_min_levels()
levels = copy.deepcopy(potential_mrca_graph.levels)
vertex_to_level_map = copy.deepcopy(potential_mrca_graph.vertex_to_level_map)
# List of non-founder individuals
non_founder_vertices = list({x // 2 for x in potential_mrca_graph.parents_map if
                             len(potential_mrca_graph.parents_map[x]) > 0})
coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filename=tree_path,
                                                               max_parent_number=2 ** 10)
coalescent_tree.remove_unary_nodes()
valid_logger = MatcherLogger(logs_directory_name=log_directory)
valid_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph, logger=valid_logger)
valid_alignments_dict = valid_matcher.find_mapping()
number_of_clades = len(valid_alignments_dict)
print(f"Clades: {number_of_clades}")
root_vertex = next(iter(valid_alignments_dict))
# Find the root vertex's spouses and form a list containing them and the root vertex
root_vertex_children_all_ploids = [
    child_ploid
    for child in potential_mrca_graph.children_map[root_vertex]
    for child_ploid in [2 * (child // 2), 2 * (child // 2) + 1]
]

correct_root_assignments = {p // 2 for child_ploid in root_vertex_children_all_ploids
                            for p in potential_mrca_graph.parents_map[child_ploid]}

os.chdir(log_directory)
save_alignment_result_to_files(valid_alignments_dict, coalescent_tree, potential_mrca_graph, tree_filename)
os.chdir("..")
valid_alignments = [d for lst in valid_alignments_dict.values() for d in lst]
valid_alignments_dict = None
# TODO: Support multiple clades
initial_average_distance_to_identity = None
identity_solution = None
initial_root_vertex_individual_assignments = None
initial_root_vertex_individual_assignments_number = 0
total_number_of_root_vertex_individual_assignments = 0
new_individual_assignments_number = 0
new_individual_simulations_number = 0

if number_of_clades == 1:
    # Supporting multiple clades requires a bit more complicated logic, so it's skipped for now
    # Calculating the average distance to identity
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
    general_result_file.write(f"Clade root individual: {root_vertex // 2}\n")
    general_result_file.write(f"Clade root and their spouses: {correct_root_assignments}\n")
    general_result_file.write(f"Individual assignments for the root: {initial_root_vertex_individual_assignments}\n")
    general_result_file.flush()

try:
    for i in range(simulation_steps):
        os.makedirs(f"{i}")
        os.chdir(f"{i}")
        # Selecting random ploids and transforming their ids into individual ids
        # random_individuals = random_subselect(non_founder_vertices, error_rate)
        random_individuals = random_subselect_poisson(non_founder_vertices, error_rate)
        print(f"Selected {len(random_individuals)} individuals")
        remove_edges = []
        add_edges = []
        for vertex in random_individuals:
            vertex_ploids = [2 * vertex, 2 * vertex + 1]
            vertex_parents = {y for x in vertex_ploids for y in potential_mrca_graph.parents_map[x]}
            vertex_parents_individuals = list({y // 2 for y in vertex_parents})
            # assert 0 < len(vertex_parents_individuals) < 3
            # assert len(vertex_parents) == 2 * len(vertex_parents_individuals)
            # vertex_parents = potential_mrca_graph.parents_map[vertex]
            random_parent = random.sample(vertex_parents_individuals, 1)[0]
            random_parent_ploids = [2 * random_parent, 2 * random_parent + 1]
            # The levels must be the same for both ploids
            random_parent_level = potential_mrca_graph.vertex_to_level_map[2 * random_parent]
            # Find the ploid id of the vertex that is connected to this parent
            random_parent_children_ploids = [x for x in vertex_ploids
                                             if potential_mrca_graph.parents_map[x] == random_parent_ploids]
            # Sometimes the same parent is specified twice as both the mother and the father
            # assert len(random_parent_children_ploids) == 1
            child_ploid = random_parent_children_ploids[0]
            # Select the new parent, start by the selecting the ploid
            new_parent_ploid_candidates = [x for x in potential_mrca_graph.levels[random_parent_level] if x
                                           not in vertex_parents]
            if not new_parent_ploid_candidates:
                warnings.warn("No other vertices at the level, skipping the error")
                continue
            new_parent_ploid = random.sample(new_parent_ploid_candidates, 1)[0]
            # new_parent_ploid = [x for x in random.sample(potential_mrca_graph.levels[random_parent_level], 5)
            #                     if x not in vertex_parents][0]
            new_parent = new_parent_ploid // 2
            # Reconnect the vertices
            # print(f"Removing {vertex}-{random_parent}, adding {vertex}-{new_parent}")
            potential_mrca_graph.remove_edge(parent=2 * random_parent, child=child_ploid, recalculate_levels=False)
            potential_mrca_graph.remove_edge(parent=2 * random_parent + 1, child=child_ploid, recalculate_levels=False)
            potential_mrca_graph.add_edge(parent=2 * new_parent, child=child_ploid, recalculate_levels=False)
            # We should recalculate levels here to be safe, but this takes a lot of time for large pedigrees!
            potential_mrca_graph.add_edge(parent=2 * new_parent + 1, child=child_ploid, recalculate_levels=False)
            add_edges.extend([(2 * random_parent, child_ploid), (2 * random_parent + 1, child_ploid)])
            remove_edges.extend([(2 * new_parent, child_ploid), (2 * new_parent + 1, child_ploid)])
        print(f"Recalculating the levels and the ancestor map")
        potential_mrca_graph.initialize_vertex_to_level_map()
        potential_mrca_graph.initialize_potential_mrca_map()

        logger = MatcherLogger(logs_directory_name=f"log_{i}")
        graph_matcher = GraphMather(coalescent_tree=coalescent_tree, processed_graph=potential_mrca_graph,
                                    logger=logger)
        result = graph_matcher.find_mapping()
        # Find the simulation solutions
        # assert len(result) == number_of_clades
        save_alignment_result_to_files(result, coalescent_tree, potential_mrca_graph, tree_filename)
        resulting_alignments = [d for lst in result.values() for d in lst]
        # incorrect_alignments = [x for x in resulting_alignments if x not in valid_alignments]
        # total_number_of_incorrect_solutions += len(incorrect_alignments)
        log_solutions(resulting_alignments, i)
        general_result_file.write("#####################################\n")
        general_result_file.write(f"Simulation step: {i}\n")
        general_result_file.write(f"Number of alignments in this simulation {len(resulting_alignments)}\n")
        # general_result_file.write(f"Total number of incorrect alignments so far among all the simulations is:"
        #                           f" {total_number_of_incorrect_solutions}\n")
        if number_of_clades == 1 and len(resulting_alignments) > 0:
            # Calculate the distances
            total_number_of_solutions += len(resulting_alignments)
            total_distance_to_solution += sum(
                get_alignments_proband_distance_probands_ignored(alignment, identity_solution, coalescent_tree.probands)
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
            total_number_of_root_vertex_individual_assignments += len(simulation_root_vertex_individual_assignments)
            new_individual_assignments_number += len({x for x in simulation_root_vertex_individual_assignments
                                                      if x not in initial_root_vertex_individual_assignments})
            if new_individual_assignments_number > 0:
                new_individual_simulations_number += 1
            # Update the counters
            if simulation_root_vertex_individual_assignments.issubset(correct_root_assignments):
                only_correct_root_assignments_counter += 1
            elif simulation_root_vertex_individual_assignments.intersection(correct_root_assignments):
                present_correct_root_assignments_counter += 1
            else:
                no_correct_root_assignments_counter += 1
        elif len(resulting_alignments) == 0:
            no_solutions_counter += 1
        general_result_file.write("#####################################\n")
        general_result_file.flush()
        # Reconstruct the initial graph
        for parent, child in add_edges:
            potential_mrca_graph.add_edge(parent=parent, child=child, recalculate_levels=False)
        for parent, child in remove_edges:
            potential_mrca_graph.remove_edge(parent=parent, child=child, recalculate_levels=False)
        potential_mrca_graph.initialize_vertex_to_level_map()
        # assert potential_mrca_graph.vertex_to_level_map == vertex_to_level_map
        # The parents and children maps are defaultdict, the dictionaries can hold the same data, but be different
        # assert not [x for x in set(initial_parents_map).union(potential_mrca_graph.parents_map)
        #             if sorted(potential_mrca_graph.parents_map[x])
        #             != sorted(initial_parents_map[x])]
        # assert not [x for x in set(initial_children_map).union(potential_mrca_graph.children_map)
        #             if sorted(potential_mrca_graph.children_map[x])
        #             != sorted(initial_children_map[x])]
        os.chdir("..")
except KeyboardInterrupt:
    print("Stop the simulation, log the final results")
# Logging the overall results
if number_of_clades == 1:
    simulations_with_solutions = simulation_steps - no_solutions_counter
    alignments_to_identity_average_distance = (total_distance_to_solution / total_number_of_solutions)
    delta = alignments_to_identity_average_distance - initial_average_distance_to_identity
    number_of_probands = len(coalescent_tree.probands)
    non_proband_vertices = len(coalescent_tree.vertex_to_level_map) - number_of_probands
    general_result_file.write(f"Initial average distance to the identity is: {initial_average_distance_to_identity}\n")
    general_result_file.write(f"Average distance from simulation solutions to the identity "
                              f"is: {alignments_to_identity_average_distance}\n")
    general_result_file.write(f"Average distance delta is: {delta}\n")
    general_result_file.write(f"The number of probands is {number_of_probands}\n")
    general_result_file.write(f"The number of non-proband vertices is: {non_proband_vertices}\n")
    general_result_file.write(f"Average distance is {alignments_to_identity_average_distance / non_proband_vertices} "
                              f"from the number of non-proband vertices in the coalescent tree\n")
    general_result_file.write("----------------------------------------------------------------------\n")
    general_result_file.write(f"The number of simulations with solutions: {simulations_with_solutions}\n")
    general_result_file.write(f"The number of simulations with no solutions: {no_solutions_counter}\n")
    general_result_file.write(f"The number of simulations with only correct assignments: "
                              f"{only_correct_root_assignments_counter}\n")
    general_result_file.write(f"The number of simulations containing correct assignments: "
                              f"{present_correct_root_assignments_counter}\n")
    general_result_file.write(f"The number of simulations not containing correct assignments: "
                              f"{no_correct_root_assignments_counter}\n")
    general_result_file.write("----------------------------------------------------------------------\n")
    general_result_file.write(f"Initial number of individual candidates for the root: "
                              f"{initial_root_vertex_individual_assignments_number}\n")
    general_result_file.write(f"Average number of individual candidates for the root per simulation with solutions: "
                              f"{total_number_of_root_vertex_individual_assignments / simulations_with_solutions}\n")
    general_result_file.write(f"New individual root assignments: {new_individual_assignments_number}\n")
    general_result_file.write(f"Average number of new individual root assignments per simulation with solutions: "
                              f"{new_individual_assignments_number / simulations_with_solutions}\n")
    general_result_file.write(f"Number of new individual root assignments simulations: "
                              f"{new_individual_simulations_number}\n")
