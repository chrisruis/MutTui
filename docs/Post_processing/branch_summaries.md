# Generating summaries from MutTui output

### Summarising an individual run

To summarise the number of mutations identified in a run, both in total and for each applied label, and to calculate the number of times each mutation has been identified, run:
```
python3 run_summary.py
```

This should be run from the directory containing the MutTui output files

Several summaries are printed to the screen. This also produces an output file all_mutation_counts.csv which contains the number of times each mutation was inferred to occur within a run. This may be useful to identify potentially problematic mutations that can be filtered out (see Mutation filtering section)

### Summarising the mutations on each branch in the tree

A set of summaries can be generated for each branch in the tree using:
```
python3 summarise_branch_mutations.py -m all_included_mutations.csv -t annotated_tree.nexus -o out_prefix
```

The all_included_mutations.csv and annotated_tree.nexus files are in the MutTui output directory

This outputs 2 files:
1) A csv file containing one row per branch with the total mutations, number of each mutation type and proportion of each mutation type on each branch
2) A nexus tree in which each branch is labelled with the total mutations, number of each mutation type and proportion of each mutation type
