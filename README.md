# CS466-Genome-Alignment-Project

One of the main uses of the sequence alignment technique is to construct phylogenetic trees,
which depicts the evolutionary relationships between different species. The amount two
sequences differ is related to the evolutionary distance between two species. This is because a
new species is created through genetic mutations, and so two species with a more common
ancestor will have less differences in their genome. In our project, we will compare the result of
an alignment between the genome sequences of a cat and a tiger and the alignment between the
genome sequences of a dog and a tiger. We expect the cat and the tiger alignment to have a
higher score and be better aligned, as a cat and tiger likely have a more recent common ancestor
than a dog and a tiger and are also in a common family (Felidae). To do this, we will implement
the dynamic programming global alignment strategy and optimize it to work on larger datasets.

Todo:
dp recursive -> memo iter
change backtracking