.. Sphinx documentation can be found here: https://www.sphinx-doc.org/en/master/index.html

================================================================================
Pedigree-ARG matching
================================================================================

This project addresses the alignment problem between a pedigree and an Ancestral Recombination Graph (ARG) that describes genetic transmissions within a given pedigree.

=====================================
Problem definition
=====================================

Specifically, given the following inputs:

1. **The Pedigree** :math:`P`;
2. **The ARG** :math:`A`: An ARG with a set of leaf vertices :math:`L`, representing genetic transmissions within :math:`P`;
3. **Initial Assignments**: A mapping :math:`f: L \to 2^{V(P)} \setminus \{\emptyset\}`, specifying initial relationships between ARG leaves and pedigree vertices,

the algorithm computes all possible **alignments** between :math:`A` and :math:`P`. An alignment is a function that assigns every vertex in the ARG to a vertex in the pedigree.

Formally, it is defined as :math:`h: V(A) \to 2^{V(P)} \setminus \{\emptyset\}`.

In other words, the goal is to find valid extensions of the initial assignments provided as input.


================================
Installation
================================

To install:

1. Clone this repository.
2. In your local copy, open a terminal.
3. Run ``pip install .``.

================================
Running the alignment
================================

The most convenient way of running this software is by using the driver file.

--------------------------------
Driver file
--------------------------------

The driver file is a YAML file that specifies the input data in the following user-friendly format:

.. code-block:: yaml

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
   output_path: "example_output"
   alignment_mode: "default"


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Initial Assignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``initial_assignments`` field specifies the initial mapping for the leaf vertices in the coalescent tree. Each object in the list must include two fields: ``coalescent_id`` and ``pedigree_ids``.

- **coalescent_id**: Specifies the ID of the coalescent node in the tree.
- **pedigree_ids**: A list of values denoting the ploids in the pedigree to which the specified ``coalescent_id`` can be mapped.

Each ploid is represented by the pedigree ID followed by one of the two ploid types:
  - **P**: Represents the paternal ploid.
  - **M**: Represents the maternal ploid.

.. note::
   The mapping can be specified to a subset of the leaf vertices. In this case, the algorithm will take the ascending
   tree for those leaf vertices.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Graph Parsing Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``pedigree`` and ``coalescent_tree`` sections describe the parsing rules used by the program. The general structure is as follows:

.. code-block:: yaml

   graph_parsing_rules:
     path: "example_coalescent_tree"
     missing_parent_notation: "-1"        # (optional)
     separation_symbol: " "               # (optional)
     skip_first_line: false               # (optional)

- **path**: Specifies the path to the file to be parsed.
  The program resolves this path using the following logic:

  1. If the path is an absolute path or a valid relative path from the current working directory (CWD), this file is used directly.
  2. Otherwise, the program treats the path as relative to the driver file's location.

Once the file is located, the program checks for any user-specified parsing rules. Below are the customizable options:

- **missing_parent_notation** *(optional)*: The sequence of characters used to denote an unknown parent. The default value is ``"-1"``.

- **separation_symbol** *(optional)*: The sequence of characters used to separate the columns in the input file. The default value is ``" "``.

- **skip_first_line** *(optional)*: Indicates whether the first line in the file should be skipped. This can be useful when the file includes a header. The default value is ``false``.

""""""""""""""""""""""""""""""
Pedigree File Structure
""""""""""""""""""""""""""""""

When parsing the pedigree file, the program assumes the following column structure:

- **First column** → Child ID
- **Second column** → Father ID
- **Third column** → Mother ID

Any additional columns in the file are ignored.

---------------------------------------------------------------

"""""""""""""""""""""""""""""""
Example Input File
"""""""""""""""""""""""""""""""

.. code-block:: text

   # id parent0 parent1
   0;;
   1;;2
   2;3;;
   3;4;5

To parse this file, use the following YAML configuration:

.. code-block:: yaml

   pedigree:
     path: "example_pedigree_path"
     missing_parent_notation: ""
     separation_symbol: ";"
     skip_first_line: true


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Alignment mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can also optionally specify the **alignment mode** that the algorithm will use by using the ``alignment_mode``
option in the input file. There are two possible modes that you can use:

1. ``default`` — Uses the default alignment approach and finds all the valid alignments between
the pedigree and the tree. This is the default mode.

2. ``example_per_root_assignment`` — Finds only one alignment per a valid ploid pedigree candidate. In other words,
it only finds a subset of all the valid alignments. This can be useful when you care only about possible clade root
assignments and want to avoid generating a large number of alignments.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``output_path`` field specifies where the results of the alignment should be stored. You can either specify an
absolute path or a relative path. If you specify a relative path, it will be interpreted relative to the
driver file's location.

---------------------------------------------------------------

As you can see, the driver file allows us to specify all the required input for the algorithm, as specified in
the Problem statement, in a customizable and a user-friendly way.

--------------------------------
Running the alignment script
--------------------------------

You can now use the driver file that you've created by running the ``scripts/driver_file/run_driver_file.py``. You can
either run the script without any arguments and get prompted by the program, or use the ``-f`` flag for specifying
the driver file's location.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The repository contains one small example for you to get started under ``scripts/driver_file/example`` which contains
a small pedigree, a coalescent tree and the driver file. The input graphs are depicted below:

.. figure:: _static/pedigree.svg
   :align: center
   :alt: Pedigree

   Example pedigree

--------------------------------------------

.. figure:: _static/tree.svg
   :align: center
   :height: 400px
   :alt: Coalescent tree

   Example coalescent tree

The initial mapping is given as follows:

- ``1`` can be mapped to any ploid of ``10`` and ``11``.
- ``2`` can be mapped to any ploid of ``12``.
- ``3`` can be mapped to any ploid of ``15``.

After running the script, you should get 8 valid alignments in the end. You can run the example by installing the
library and running the following command from the root directory:

.. code-block:: bash

   python scripts/driver_file/run_driver_file.py -f "scripts/driver_file/example/driver_file.yaml"


--------------------------------
Output format
--------------------------------

The algorithm saves the alignment results in the directory specified by ``output_path``.
For every clade in the coalescent tree, a separate folder will be created, named after the clade root's ID.

Inside each clade folder:

1. **Alignment Files**:

   - The folder will contain a set of alignment files, where each file represents a distinct alignment.

   - Two alignments are considered distinct if there is at least one coalescent vertex that is mapped differently.

   - Each alignment file will be named as ``alignment_{id}``, where ``{id}`` is a sequential counter identifying the alignment.

2. **Statistics File**:

   - A file named ``clade_{root_id}`` will also be included. This file aggregates the data and provides statistical insights about the alignments.

.. note:: **No Alignments**

    If no alignments are generated, it indicates that the input data contains errors, and no solutions exist for the given input.

.. note:: **Alignment Similarity**

    A high number of alignments does not necessarily mean high variability. Many alignments might be very similar.
    To evaluate the true diversity of alignments, refer to the statistics file, which provides the number of
    possibilities for each coalescent vertex.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Alignment Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each alignment result is stored as a separate plain text file with a simple structure.
The file consists of two sections:

1. **Statistical Data** — A block containing relevant statistics.
2. **Alignment Data** — A mapping of coalescent vertex IDs to pedigree ploid IDs.

The alignment data follows the format ``{coalescent_vertex_id}: {pedigree_individual_id}{ploid_type}``

- **coalescent_vertex_id** — The coalescent vertex ID.
- **pedigree_individual_id** — The individual's ID to which the coalescent vertex is mapped.
- **ploid_type** — The ploid type.

**Example File Structure:**

.. code-block:: text

    // Statistical data
    ...
    // Alignment
    5: 1M
    3: 15M
    4: 3P
    2: 12M
    1: 10P



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   simulating_errors
   error_alignments
   api