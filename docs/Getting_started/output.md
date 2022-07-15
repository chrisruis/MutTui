# MutTui output files

#### mutational_spectrum_label_X.csv

The inferred SBS spectrum from branches with label X. If your tree wasn't labelled, there will be one output file with label A. If your tree was labelled, there will be one output file per label

Column 1 contains each contextual mutation (named with **upstream base [ancestral base > mutated base] downstream base** so A[C>T]G is the number of C>T mutations with A upstream and G downstream)

Column 2 contains the count of the contextual mutation

#### mutational_spectrum_label_X.pdf

Plot of the SBS spectrum in mutational_spectrum_label_X.csv

#### DBS_label_X.csv

The inferred DBS spectrum from branches with label X. If your tree wasn't labelled, there will be one output file with label A. If your tree was labelled, there will be one output file per label

Column 1 contains each double substitution (named with **ancestral_doublet > mutated_doublet** so AC>CG is the number of mutations from an AC doublet to a CG doublet)

Column 2 contains the count of the double substitution

#### DBS_label_X.pdf

Plot of the DBS spectrum in DBS_label_X.csv

#### mutation_types_label_A.csv

The inferred spectrum of mutation counts without their context. If your tree wasn't labelled, there will be one output file with label A. If your tree was labelled, there will be one output file per label

Column 1 contains the mutation type (named with ancestral base followed by mutated base so CA is the number of C>A mutations)

Column 2 contains the count of the mutation type

#### mutation_types_label_A.pdf

Plot of the mutation type spectrum in mutation_types_label_A.csv

#### all_included_mutations.csv

All of the mutations in any of the output SBS spectra. The columns are:
* Mutation_in_alignment - the mutation and its position in the alignment
* Mutation_in_genome - the mutation and its position in the genome
* Substitution - the contextual mutation, written as "upstream base [ancestral base > mutated base] downstream base"
* Branch - the branch on which the mutation occurred. You can identify the branch in the tree by loading annotated_tree.nexus into FigTree

#### all_included_double_substitutions.csv

All of the mutations in any of the output DBS spectra. The columns are:
* Mutation_in_alignment - the mutation and its position in the alignment
* Mutation_in_genome - the mutation and its position in the genome
* Substitution - the double substitution, written as "ancestral doublet > mutated doublet"
* Branch - the branch on which the mutation occurred. You can identify the branch in the tree by loading annotated_tree.nexus into FigTree
* Original_mutation - specifies whether the original mutation was on the forward or reverse strand

#### mutations_not_included.csv

Mutations inferred by treetime that were not included in an output SBS or DBS spectrum and the reason they were excluded

#### gaps_to_N_alignment.fasta

The input alignment with gaps converted to Ns. This is used because gaps are treated as characters and included in reconstructions by treetime while Ns are treated as missing data

# treetime output files

These files are included in the MutTui output directory but are produced by treetime during ancestral reconstruction

#### branch_mutations.txt

Contains the mutations inferred by treetime through ancestral reconstruction

#### ancestral_sequences.fasta

Alignment from treetime containing observed and inferred ancestral sequences

#### annotated_tree.nexus

Tree from treetime containing internal node labels and inferred mutations on each branch

#### sequence_evolution_model.txt

Inferred model parameters from treetime ancestral reconstruction.
