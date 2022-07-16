# Splitting a single spectrum into multiple groups

The output spectrum from MutTui can be split post-hoc into multiple groups using a branch labels file:
```
python3 post_process_branch_mutations.py -m all_included_mutations.csv -l branch_labels.csv -o split_out
```

all_included_mutations.csv is in the MutTui output directory and contains all of the mutations included in an output spectrum and their branch. The branch_labels.csv file contains branch names and their label. MutTui outputs one spectrum per label in branch_labels.csv that contains all of the mutations on the branches with that label
