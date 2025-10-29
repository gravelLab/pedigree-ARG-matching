# Pedigree-ARG Matching

This project addresses the alignment problem between a pedigree and a coalescent tree that describes genetic transmissions within a given pedigree. Specifically, given the following inputs:

1. **The Pedigree $P$**;
2. **The coalescent tree $T$**: A coalescent tree with a set of leaf vertices $L$, representing genetic transmissions within $P$;
3. **Initial Assignments**: A mapping $`f: L \to 2^{V(P)} \setminus \{\emptyset\}`$, specifying initial relationships between ARG leaves and pedigree vertices,

the algorithm computes all possible **vertex alignments** between $T$ and $P$. 

A **vertex alignment** is a function that assigns each vertex in $T$ to a vertex in $P$. That is:

$$
h: V(T) \to V(P)
$$

We say that a vertex alignment is **valid** if the corresponding vertices in $P$ represent a history
of genetic transmissions that respect $T$.

In order to verify that a given vertex alignment is indeed valid,
the algorithms finds at least one valid **edge alignment** (that is, an alignment that maps
every edge in $E(T)$ to a path in $P$) that corresponds to the given vertex alignment.

In other words, the goal is to find valid extensions of the initial assignments provided as input.

## Documentation

The documentation for the project can be found [here](https://gravellab.github.io/pedigree-ARG-matching/).

## Installation

1. Clone the repository and `cd` into the root of the project.
2. Run `pip install .`
