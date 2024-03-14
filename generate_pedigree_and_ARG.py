import os

import msprime.pedigrees

import msprime.ancestry
from tskit import TreeSequence, TableCollection, Tree
from genealogical_graph import CoalescentTree

sequence_length = 2
recombination_rate = 0.1
arg_filename = "arg.arg"


def generate_and_save_ARG(parsed_pedigree: TableCollection):
    simulated_arg: TreeSequence = msprime.ancestry.sim_ancestry(initial_state=parsed_pedigree,
                                                                model="fixed_pedigree",
                                                                recombination_rate=recombination_rate,
                                                                sequence_length=sequence_length)
    coalescent_tree_filename = arg_filename
    coalescent_tree: CoalescentTree = CoalescentTree.get_arg(tree_sequence=simulated_arg)
    coalescent_tree.save_to_file(filename=coalescent_tree_filename)
    simulated_arg.draw_svg(path=f"{coalescent_tree_filename}.svg", size=(1000, 1000))
    return simulated_arg


def generate_and_save_coalescent_tree(tree_counter: int, tree: Tree):
    coalescent_tree_filename = f"coalescent_tree_{tree_counter}"
    tree_counter += 1
    coalescent_tree = CoalescentTree.get_coalescent_tree(tree)
    coalescent_tree.save_to_file(filename=coalescent_tree_filename)
    tree.draw_svg(path=f"{coalescent_tree_filename}.svg", size=(1000, 1000))
    return tree_counter


population_size = 1004
end_time = 5

# Create a directory for the files
directory_path = f"pedigrees/{population_size}_{end_time}"

if not os.path.exists(directory_path):
    os.makedirs(directory_path)
os.chdir(directory_path)

seed = 1234

filename_prefix = f"{population_size}_{end_time}"
filename = f"{filename_prefix}.pedigree"
file = open(filename, 'w')
pedigree = msprime.pedigrees.sim_pedigree(population_size=population_size,
                                          random_seed=seed,
                                          sequence_length=sequence_length,
                                          end_time=end_time)
pedigree: TableCollection
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

arg: TreeSequence = generate_and_save_ARG(pedigree)
print("Generated ARG\n")
first_tree: Tree = arg.first()

generate_and_save_coalescent_tree(0, first_tree)
arg_parsed: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filename=arg_filename,
                                                                          max_parent_number=2 ** 10)
