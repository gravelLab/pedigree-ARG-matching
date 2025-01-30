================================================================================
Simulating errors
================================================================================

The developed algorithm always finds solutions when there are no errors in the input data. However, in a more realistic
scenarios, we will most likely have errors both in the pedigree and the ARG. In order to test how stable towards
potential errors the developed algorithm is, the library contains functionality for simulating errors both in a
pedigree and a coalescent tree.

The general idea behind this approach is as follows:

1. Get perfect data for the pedigree and the coalescent tree using msprime.
2. Simulate errors either within the pedigree or the tree.
3. Run multiple alignments with the modified data and compare the alignment results.

More information on this approach can be found in :doc:`error_alignments`.

-----------------------------------------------
Simulating pedigree errors
-----------------------------------------------

Some information regarding how pedigree errors are simulated

-----------------------------------------------
Simulating coalescent tree errors
-----------------------------------------------

Some information regarding how tree errors are simulated
