# MutTui output files

### mutational_spectrum_label_X.csv

The inferred mutational spectrum along branches labelled with label X. There will be one of these files for each branch label with the "X" part of the file name replaced with the label. Note that if your tree does not contain multiple labels (i.e. you ran MutTui without -lt and -l), you will have one output mutational spectrum called mutational_spectrum_label_A.csv.

In this file, column 1 contains each contextual mutation (96 mutations with DNA, 192 with RNA). The format of these contextual mutations is:
* upstream base [ancestral base > mutated base] downstream base
So the mutation is from "ancestral base" to "mutated base". The "upstream base" and "downstream base" give the context of the mutation.

Column 2 contains the number of the contextual mutation inferred to have occured within the evolutionary history of the dataset.

### mutational_spectrum_label_X.pdf

Plot of the mutational spectrum in mutational_spectrum_label_X.csv. Each of the output spectra will have a plot.

### all_included_mutations.csv

All of the mutations that were included in any of the output mutational spectra are included in this file. The columns are:
* Mutation_in_alignment - the mutation and its position in the alignment as ancestral base followed by alignment position followed by mutated base
* Mutation_in_genome - the mutation and its position in the genome as ancestral base followed by genome position followed by mutated base
* Substitution - the contextual mutation, written as "upstream base [ancestral base > mutated base] downstream base"
* Branch - the branch the mutation occurred along. You can identify the branch in the tree by loading annotated_tree.nexus into a tree viewing programme (e.g. FigTree)
* Label - the label of the branch along which the mutation occurred

### mutations_not_included.csv

Inferred mutations that were not included in the output spectra. Mutations are not included if they are at the first or last positions in the genome (as there is no context on one side of these positions) or if three or more adjacent genome positions change along the same branch, as we cannot infer whether these mutations occurred independently, as double substitutions or as triple or greater substitutions.

### ancestral_sequences.fasta

An alignment from treetime containing each of the sequences provided to MutTui with -a and the inferred nucleotide at each internal node at each alignment position. The names of the ancestral sequences can be reconciled with the tree by comparing with annotated_tree.nexus.

### annotated_tree.nexus

Tree from treetime. This will be the same as the tree that was provided to MutTui with -t but each node in the tree is named with its node label and each branch in the tree is labelled with the mutations that occurred along the branch. You can view this tree with FigTree and view the mutations along each branch by switching on branch labels and labelling by mutations. This will be useful if you want to check when particular mutations occurred. MutTui uses this tree to identify mutations.

### sequence_evolution_model.txt

Inferred model parameters from treetime ancestral reconstruction.
