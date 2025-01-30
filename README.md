# Pedigree-ARG Matching

This project addresses the alignment problem between a pedigree and an Ancestral Recombination Graph (ARG) that describes genetic transmissions within a given pedigree. Specifically, given the following inputs:

1. **The Pedigree $P$**;
2. **The ARG $A$**: An ARG with a set of leaf vertices $L$, representing genetic transmissions within $P$;
3. **Initial Assignments**: A mapping $f: V(L) \to 2^{P} \setminus \emptyset\$, specifying initial relationships between ARG leaves and pedigree vertices,

the algorithm computes all possible **alignments** between $A$ and $P$. An alignment is a function that assigns every vertex in the ARG to a vertex in the pedigree. 
Formally, it is defined as $h: V(A) \to 2^{P} \setminus \emptyset$. 
In other words, the goal is to find valid extensions of the initial assignments provided as input.

## Documentation

The documentation for the project can be found here: 

## Installation

1. Clone the repository and `cd` into the root of the project.
2. Run `pip install .`
