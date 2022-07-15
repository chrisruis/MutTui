# Quick start

## Basic usage

To calculate a mutational spectrum, 4 files are needed:
1) Sequence alignment in fasta format
2) Phylogenetic tree in newick format
3) Reference genome in fasta format
4) Conversion from alignment position to genome position

Collect these files and then run:

```
mkdir muttui_out
python3 muttui.py -a alignment.fasta -t tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```

## Creating the position conversion file

To create the conversion.txt file from a VCF, run:

```
python3 convert_vcf_position_translation.py -v positions.vcf -o conversion.txt
```

## Labelling your tree

More details in [Labelling_your_tree](https://github.com/chrisruis/MutTui/tree/main/docs/Labelling_your_tree)

By default, MutTui will assign all branches in the given tree to the same mutational spectrum. If you want to split your tree into multiple groups, first label the tree using:

```
python3 add_tree_node_labels.py -t unlabelled_tree.nwk -o tree_with_node_labels.nwk
python3 label_tree.py -t unlabelled_tree.nwk -r root_state -s branch_state branch_state -o state_labelled_tree.nwk
```

Then run MutTui including the labelled tree with -lt:

```
mkdir muttui_out
python3 muttui.py -a alignment.fasta -t tree.nwk -lt labelled_tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```
