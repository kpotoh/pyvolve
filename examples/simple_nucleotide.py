# This example script demonstrates how to evolve according to a simple nucleotide model. All model parameters are default: equal mutation rates and equal equilibrium frequencies (e.g. JC69 model).

import pyvolve

# Define a phylogeny, from a file containing a newick tree
my_tree = pyvolve.read_tree(file = "file_with_tree.tre")

# Define a nucleotide model, as a Model object.
my_model = Model("nucleotide")

# Assign the model to a Partition. The size argument indicates to evolve 250 positions
my_partition = Partition(models = my_model, size = 250)

# Evolve!
my_evolver = Evolver(partitions = my_partition, tree = my_tree)
my_evolver()