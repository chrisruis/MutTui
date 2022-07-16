# Generating branch summaries from MutTui output

A set of summaries can be generated for each branch in the tree using:
```
python3 summarise_branch_mutations.py -m all_included_mutations.csv -t annotated_tree.nexus -o out_prefix
```

The all_included_mutations.csv and annotated_tree.nexus files are in the MutTui output directory

This outputs 2 files:
1) A csv file containing one row per branch with the total mutations, number of each mutation type and proportion of each mutation type on each branch
2) A nexus tree in which each branch is labelled with the total mutations, number of each mutation type and proportion of each mutation type
