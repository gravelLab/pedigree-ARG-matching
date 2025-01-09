# Pedigree-ARG Matching

This project addresses the alignment problem between a pedigree and an Ancestral Recombination Graph (ARG) that describes genetic transmissions within a given pedigree. Specifically, given the following inputs:

1. **The Pedigree $P$**;
2. **The ARG $A$**: An ARG with a set of leaf vertices $L$, representing genetic transmissions within $P$;
3. **Initial Assignments**: A mapping $f: V(L) \to 2^{P} \setminus \emptyset\$, specifying initial relationships between ARG leaves and pedigree vertices,

the algorithm computes all possible **alignments** between $A$ and $P$. An alignment is a function that assigns every vertex in the ARG to a vertex in the pedigree. 
Formally, it is defined as $h: V(A) \to 2^{P} \setminus \emptyset$. 
In other words, the goal is to find valid extensions of the initial assignments provided as input.

## Initial Alignment

The algorithm requires the **initial alignment** to be specified in a YAML format as follows:

```yaml
initial_assignments:
  - coalescent_id: 1
    pedigree_ids: [11P]
  - coalescent_id: 2
    pedigree_ids: [11M]
  - coalescent_id: 3
    pedigree_ids: [12P, 13M, 13P]
  - coalescent_id: 4
    pedigree_ids: [14M, 15P, 19M]
coalescent_tree_path: "path/to/coalescent/tree"
pedigree_path: "path/to/pedigree/file"
```

Here, every object of the ```initial_assignments``` list specifies the list of pedigree values for the specified 
```coalescent_id```. Every pedigree id is specified in the format ```{individual_id}{ploid_type}``` where 
```individual_id``` is the id of the pedigree individual (which must be a valid id in the pedigree) and the ploid type 
is either "M" or "P" where "M" specifies the maternal ploid () and "P" specifies the paternal ploid.
The values for ```pedigree_path``` and ```coalescent_tree_path``` specify the paths to the pedigree and the coalescent 
tree  respectively.


