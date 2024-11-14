import warnings

from alignment.graph_matcher import *
from scripts.utility import *

simulation_root_dir = "error_pedigrees"

error_rate = 0.01

os.makedirs(simulation_root_dir, exist_ok=True)
os.chdir(simulation_root_dir)

simulation_name = get_non_existing_directory_path("Specify the name of the simulation: ")
simulation_steps = get_natural_number_input("Specify the number of simulation steps: ")
pedigree_path = get_file_path("Specify the absolute path to the pedigree:")

os.makedirs(simulation_name)
os.chdir(simulation_name)

potential_mrca_graph = PotentialMrcaProcessedGraph.get_processed_graph_from_file(filepath=pedigree_path,
                                                                                 initialize_levels=True,
                                                                                 initialize_ancestor_maps=False)
# List of non-founder individuals
non_founder_vertices = list({x // 2 for x in potential_mrca_graph.parents_map})

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
        for index, vertex in enumerate(random_individuals):
            print(f"Processing {index} / {len(random_individuals)}")
            vertex_ploids = {2 * vertex, 2 * vertex + 1}
            vertex_parents = {y for x in vertex_ploids for y in potential_mrca_graph.parents_map.get(x, [])}
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
                                             if potential_mrca_graph.parents_map.get(x, []) == random_parent_ploids]
            # Sometimes the same parent is specified twice as both the mother and the father
            # assert len(random_parent_children_ploids) == 1
            child_ploid = random_parent_children_ploids[0]
            # Select the new parent, start by the selecting the ploid
            max_attempts = len(potential_mrca_graph.levels[random_parent_level]) * 2
            attempts = 0
            while attempts < max_attempts:
                new_parent_ploid = random.choice(potential_mrca_graph.levels[random_parent_level])
                if new_parent_ploid not in vertex_parents:
                    break
                attempts += 1
            else:
                new_parent_ploid_candidates = [x for x in potential_mrca_graph.levels[random_parent_level] if x
                                               not in vertex_parents]
                if not new_parent_ploid_candidates:
                    warnings.warn("No other vertices at the level, skipping the error")
                    continue
                new_parent_ploid = random.sample(new_parent_ploid_candidates, 1)[0]
            new_parent = new_parent_ploid // 2
            # Reconnect the vertices
            # print(f"Removing {vertex}-{random_parent}, adding {vertex}-{new_parent}")
            potential_mrca_graph.remove_edge(parent=2 * random_parent, child=child_ploid, recalculate_levels=False)
            potential_mrca_graph.remove_edge(parent=2 * random_parent + 1, child=child_ploid, recalculate_levels=False)
            potential_mrca_graph.add_edge(parent=2 * new_parent, child=child_ploid, recalculate_levels=False)
            potential_mrca_graph.add_edge(parent=2 * new_parent + 1, child=child_ploid, recalculate_levels=False)
            add_edges.extend([(2 * random_parent, child_ploid), (2 * random_parent + 1, child_ploid)])
            remove_edges.extend([(2 * new_parent, child_ploid), (2 * new_parent + 1, child_ploid)])
        potential_mrca_graph.save_to_file(filepath=f"{simulation_name}_{i}.pedigree")
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
