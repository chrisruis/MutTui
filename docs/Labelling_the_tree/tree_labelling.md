# Why label your tree?

It's often useful to compare spectra between different groups of pathogens. If these groups are separate clades, MutTui can be run on each clade separately and the outputs compared. But it's often useful to compare multiple groups within the same tree. For example, we might want to compare the green and purple groups within this tree:

<img src="_figures/labelled_tree.png" width = "500">

In this tree, the purple sequences form a single cluster. But it might be that the groups are more interspersed:

<img src="_figures/labelled_tree2.png" width = "500">

Or there might be several groups that are very interspersed across the tree:

<img src="_figures/labelled_tree_complex.png" width = "500">

We can calculate the spectra of these groups separately in any of these cases by labelling the tree

# How to label the tree

We label the tree in 2 steps. First add names to all branches in the tree using:
```
python3 add_tree_node_labels.py -t tree.nwk -o tree_node_names.nwk
```

This outputs the same tree given with -t but each of the internal nodes now has a label "Node" followed by a number. We need to identify the node labels on which we want the label to change. To do this, open tree_node_names.nwk in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), switch node labels on and change display to label

Once you've identified the branches on which the label should change, label the tree using:
```
python3 label_tree.py -t tree.nwk -r R -s Node10____S1 Node50____S2 Node75____S1 TipSequence1____S3 -o labelled_tree.nwk
```

This takes the same input tree as before. Option -r specifies the label at the root of the tree. With -s, we specify each of the branches on which we want the label to change and what we want the label to change to. Separate the branch from the label with four underscores. We can specify as many label changes as we want. Here, we've specified 4 label changes and, for example, has specified that the label should change to "S1" along the branch named Node10. To change the label on a tip branch, give the name of the tip sequence

If there are a large number of label changes, rather than specifying them all with -s, you can provide a file with -sf that contains the branches on which the label changes and the label to change to

You can check that the labelling has worked as expected by opening labelled_tree.nwk in FigTree and searching for four underscores followed by the label (to check the tip labels) and by switching node labels on (to check the internal node labels)

The output of label_tree.py is a labelled tree which can be provided to MutTui with -lt, for example:
```
python3 muttui.py -a alignment.fasta -t tree.nwk -lt labelled_tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```

This will output one spectrum per label

# Inferring state changes as part of the MutTui run

It's usually best to identify where the label should change and label the tree before running MutTui as above. However, there may be cases where this isn't possible but the label of each of the tips is known. In these cases, you can infer where the label changes across the tree as part of the MutTui run using:

```
python3 muttui.py -a alignment.fasta -t tree.nwk -l tip_states.csv -r reference.fasta -c conversion.txt -o muttui_out
```

-l contains tip sequences and their labels. Any tips not included in this file are given the label OTHER (or OTHER_A if OTHER is already used in the -l file labels). MutTui then runs the mugration model within treetime to reconstruct the label onto the tree and splits branches based on this reconstruction. Again, one spectrum is output per label
