import os

import msprime.pedigrees

import msprime.ancestry
from tskit import TreeSequence
from genealogical_graph import CoalescentTree

sequence_length = 4


def generate_coalescent_tree(parsed_pedigree, tree_counter):
    # parsed_pedigree = msprime.parse_pedigree(text_file=pedigree_file, sequence_length=20).tree_sequence()

    simulated_coalescent_tree = msprime.ancestry.sim_ancestry(initial_state=parsed_pedigree,
                                                              model="fixed_pedigree",
                                                              sequence_length=sequence_length)
    simulated_coalescent_tree: TreeSequence
    for tree in simulated_coalescent_tree.trees():
        coalescent_tree_file = open(f"coalescent_tree_{tree_counter}", 'w')
        tree_counter += 1
        coalescent_tree = CoalescentTree.get_coalescent_tree(tree)
        coalescent_tree: CoalescentTree
        for level in coalescent_tree.levels:
            for vertex in level:
                if vertex in coalescent_tree.parents_map:
                    [parent] = coalescent_tree.parents_map[vertex]
                    coalescent_tree_file.write(f"{vertex} {parent}\n")
        coalescent_tree_file.close()
    return tree_counter


population_size = 1001
end_time = 16

# Create a directory for the files
directory_path = f"pedigrees/{population_size}_{end_time}"

if not os.path.exists(directory_path):
    os.makedirs(directory_path)
os.chdir(directory_path)

seed = 1234

filename = f"{population_size}_{end_time}.pedigree"
file = open(filename, 'w')
pedigree = msprime.pedigrees.sim_pedigree(population_size=population_size,
                                          random_seed=seed,
                                          sequence_length=sequence_length,
                                          end_time=end_time)
file.write("# id parent0 parent1 time\n")
for i in range(len(pedigree.individuals)):
    [left_parent, right_parent] = pedigree.individuals[i].parents
    if left_parent == "-1":
        left_parent = '.'
    if right_parent == "-1":
        right_parent = '.'
    file.write(f"{i} {left_parent} {right_parent}\n")

print("Generated pedigree\n")

file.close()

file = open(filename)

counter = 0
for i in range(sequence_length):
    counter = generate_coalescent_tree(pedigree, counter)
