# Quick start

## Basic usage

To calculate the mutational spectrum of an input alignment and tree run:

```
mkdir muttui_out
python3 muttui.py -a alignment.fasta -t tree.nwk -r reference_genome.fasta -c alignment_to_genome_position.txt -o muttui_out
```

## Labelling your tree

By default, MutTui will assign all of branches in the given tree to the same mutational spectrum. If you want to split your tree into multiple groups, run these 2 scripts before running MutTui (more details in the Labelling_your_tree section):

```
python3 add_tree_node_labels.py -t unlabelled_tree.nwk -o tree_with_node_labels.nwk
python3 label_tree.py -t unlabelled_tree.nwk -r root_state -s branch_state branch_state -o state_labelled_tree.nwk
```

You can then run MutTui splitting the branches based on the state labelled tree from label_tree.py:

```
mkdir muttui_out
python3 muttui.py -a alignment.fasta -t tree.nwk -lt state_labelled_tree.nwk -r reference_genome.fasta -c alignment_to_genome_position.txt -o muttui_out
```
