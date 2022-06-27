## Why label your tree?

Mutational spectra contain the history of contextual mutations within a given dataset. A single mutational spectrum can be very informative but it's often useful to compare mutational spectra between different phylogenetic groups, for example to investigate the mutagens or repair processes that are operating differently between them

We can compare mutational spectra with MutTui in 2 ways:
1) Calculate spectra in separate MutTui runs and then compare - use this if the different groups can't be combined into a single phylogenetic tree
2) Label the tree to calculate the spectra of different groups separately in a single MutTui run - use this if the different groups are in the same phylogenetic tree

If you want to use option 1., you can run MutTui with its default options on each clade

For option 2., we first need to label the tree

## What labelling option should I use?

* If you want to calculate the mutational spectrum of a single clade without splitting into multiple spectra:

```
python3 muttui.py -a alignment.fasta -t tree.nwk -r reference.fasta -c alignment_to_genome_conversion.txt -o muttui_out
```

This will produce a single output spectrum for your clade.

* If you know where the state changes are in your tree:

```
python3 add_tree_node_labels.py -t unlabelled_tree.nwk -o tree_with_node_labels.nwk
python3 label_tree.py -t unlabelled_tree.nwk -r root_state -s branch_state branch_state -o state_labelled_tree.nwk
python3 muttui.py -a alignment.fasta -t tree.nwk -lt state_labelled_tree.nwk -r reference_genome.fasta -c alignment_to_genome_position.txt -o muttui_out
```

This will enable you to identify the branches along which your state changes, label all branches with their state and with the state changing along your chosen branch(es) and then run MutTui splitting into the multiple states.

* If you know what each of the tip states is but you're not sure where the state changes are on the tree:

```
python3 muttui.py -a alignment.fasta -t tree.nwk -l tip_states.csv -r reference_genome.fasta -c alignment_to_genome_position.txt -o muttui_out
```

This will run the mugration model in treetime to identify the branches along which the state changes, then run MutTui to calculate the mutational spectrum of each state.
