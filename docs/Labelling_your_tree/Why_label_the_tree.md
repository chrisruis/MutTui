### Why label your tree?

Mutational spectra contain the history of contextual mutations within a given dataset. A single mutational spectrum can be very informative and be used to infer some of the underlying mutagenic and repair processes that have contributed to a mutational history. In its default mode, MutTui calculates a single mutational spectrum across the entire provided tree.

Its often very useful to compare mutational spectra between different parts of a pathogen's evolutionary history. This can be used to infer the mutagenic processes and repair mechanisms acting in different groups of pathogens and therefore identify processes that are active in some groups but not in others. To compare mutational spectra between groups, you can run MutTui separately on the two or more groups and then compare the calculated spectra. This is the best strategy if you're comparing two groups of pathogens that are relatively distantly related, for example different genotypes or sequence types within a species.

However, you might want to compare closely related clades, including in situations where one clade is nested within the diversity of another. In these cases, you can obtain the mutational spectra of the multiple clades in a single MutTui run by labelling your phylogenetic tree.

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
