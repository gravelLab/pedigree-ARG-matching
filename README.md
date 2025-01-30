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
coalescent_tree:
  path: "example_coalescent_tree"
  missing_parent_notation: "-1"
  separation_symbol: " "
  skip_first_line: false
pedigree:
  path: "example_pedigree.pedigree"
  missing_parent_notation: "-1"
  separation_symbol: " "
  skip_first_line: false
output_path: "example"
```

Here, every object of the ```initial_assignments``` list specifies the list of pedigree values for the specified 
```coalescent_id```. Every pedigree id is specified in the format ```{individual_id}{ploid_type}``` where 
```individual_id``` is the id of the pedigree individual (which must be a valid id in the pedigree) and the ploid type 
is either "M" or "P" where "M" specifies the maternal ploid () and "P" specifies the paternal ploid.
The values for ```pedigree``` and ```coalescent_tree``` specify the paths and the parsing rules for the pedigree and 
the coalescent tree  respectively.

## Running

In order to simplify the process of specifying the input for the program, the user can create a separate YAML driver 
file to tell the algorithm which input data it should use.

### Driver file specification

The driver file specifies all the information required by the algorithm. 

#### Graph Parsing Rules

The `pedigree` and `coalescent_tree` sections describe the parsing rules used by the program. The general structure is as follows:

```yaml
graph_parsing_rules:
  path: "example_coalescent_tree"
  missing_parent_notation (optional): "-1"
  separation_symbol (optional): " "
  skip_first_line (optional): false
```

- **`path`**: Specifies the path to the file to be parsed. The program resolves this path using the following logic:
  1. If the path is an absolute path or a valid relative path from the current working directory (CWD), the file is used directly.
  2. Otherwise, the program treats the path as relative to the driver file's location.

Once the file is located, the program checks for any user-specified parsing rules. Below are the customizable options:

- **`missing_parent_notation`**: The sequence of characters used to denote an unknown parent. The default value is `"-1"`.
- **`separation_symbol`**: The sequence of characters used to separate the columns in the input file. The default value is `" "`.
- **`skip_first_line`**: Indicates whether the first line in the file should be skipped. This is useful when the file includes a header. The default value is `false`.

The program assumes:
- The **first column** represents the child ID.
- The **second column** represents the father ID.
- The **third column** represents the mother ID.

Any additional columns in the file will be ignored.

---

For example, consider the following input file:

```txt
# id parent0 parent1
0;;
1;;2
2;3;;
3;4;5
```

To parse this file, you would specify the following rules:

```yaml
pedigree:
  path: "example_pedigree_path"
  missing_parent_notation: ""
  separation_symbol: ";"
  skip_first_line: true
```

---

#### Initial Assignments

The `initial_assignments` field specifies the initial mapping for the leaf vertices in the coalescent tree. Each object in the list must include two fields: `coalescent_id` and `pedigree_ids`.  

- **`coalescent_id`**: Specifies the ID of the coalescent node in the tree.  
- **`pedigree_ids`**: A list of values denoting the ploids in the pedigree to which the specified `coalescent_id` can be mapped. Each ploid is represented by the pedigree ID followed by one of the two ploid types:  
  - **`P`**: Represents the paternal ploid.  
  - **`M`**: Represents the maternal ploid.  

Note: All leaf vertices, and only the leaf vertices, must have their initial mappings specified.

---

#### Output

The `output_path` field specifies where the results of the alignment should be stored.

## Output Format

The algorithm saves the alignment results in the directory specified by `output_path`. For every clade in the coalescent tree, a separate folder will be created, named after the clade root's ID.

Inside each clade folder:

1. **Alignment Files**:
   - The folder will contain a set of alignment files, where each file represents a distinct alignment.
   - Two alignments are considered distinct if there is at least one coalescent vertex that is mapped differently.
   - Each alignment file will be named as `alignment_{id}`, where `{id}` is a sequential counter identifying the alignment.

2. **Statistics File**:
   - A file named `clade_{root_id}` will also be included.
   - This file aggregates the data and provides statistical insights about the alignments.

### Important Notes:
- **No Alignments**:
  If no alignments are generated, it indicates that the input data contains errors, and no solutions exist for the given input.

- **Alignment Similarity**:
  A high number of alignments does not necessarily mean high variability. Many alignments might be very similar. To evaluate the true diversity of alignments, refer to the statistics file, which provides the number of possibilities for each coalescent vertex.
