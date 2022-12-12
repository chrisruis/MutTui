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
MutTui run -a alignment.fasta -t tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```

## Creating the position conversion file

To create the conversion.txt file from a VCF, run:

```
MutTui convert-vcf -v positions.vcf -o conversion.txt
```

## Labelling your tree

More details in [Labelling_your_tree](https://github.com/chrisruis/MutTui/tree/main/docs/Labelling_the_tree)

By default, MutTui will assign all branches in the given tree to the same mutational spectrum. If you want to split your tree into multiple groups, first label the tree using:

```
MutTui label-nodes -t unlabelled_tree.nwk -o tree_with_node_labels.nwk
MutTui label-tree -t unlabelled_tree.nwk -r root_state -s branch_state branch_state -o labelled_tree.nwk
```

Then run MutTui including the labelled tree with -lt:

```
mkdir muttui_out
MutTui run -a alignment.fasta -t tree.nwk -lt labelled_tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```
