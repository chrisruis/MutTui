# Why label your tree?

It's often useful to compare spectra between different groups of pathogens. If these groups are separate clades, MutTui can be run on each clade separately and the outputs compared. But it's often useful to compare multiple groups within the same tree. For example, we might want to compare the green and purple groups within this tree:

<img src="_figures/labelled_tree.png">

In this tree, the purple sequences form a single cluster. But it might be that the groups are more interspersed:

<img src="_figures/labelled_tree2.png">

Or there might be several labels that very interspersed across the tree:

<img src="_figures/labelled_tree_complex.png">

We can calculate the spectra of these groups separately in a single MutTui run by labelling the tree

# How to label the tree

We can label the tree in 2 steps. The first is to add names to all branches in the tree. The second is to use these names to label all branches

We add names to all branches using:
```
python3 add_tree_node_labels.py -t tree.nwk -o tree_node_names.nwk
```

This outputs the same tree given with -t but each of the internal nodes now has a label Node followed by a number. We need to identify the node labels on which we want the label to change. To do this, open tree_node_names.nwk in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), switch node labels on and change display to label.

## Inferring state changes as part of the MutTui run

If you know what each of the tip states is but you're not sure where the state changes are on the tree:

```
python3 muttui.py -a alignment.fasta -t tree.nwk -l tip_states.csv -r reference_genome.fasta -c alignment_to_genome_position.txt -o muttui_out
```

This will run the mugration model in treetime to identify the branches along which the state changes, then run MutTui to calculate the mutational spectrum of each state.
