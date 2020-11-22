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

## Combining mutational spectra

If you want to combine multiple spectra from MutTui into a single combined spectrum (more details in the Combining_spectra section), run:

```
python3 combine_spectra.py -s spectrum1.csv spectrum2.csv spectrum3.csv spectrum4.csv -o combined_spectrum.csv
```

## Comparing mutational spectra

To compare multiple mutational spectra from MutTui (more details in the Comparing_spectra section), run:

```
python3 compare_spectra.py -s spectrum1.csv spectrum2.csv -o comparison_prefix
```

## Plotting mutational spectra

MutTui will plot each spectrum that it calculates. If you want to plot a spectrum separately (e.g. a combined spectrum) use:

```
python3 plot_spectrum.py -s spectrum.csv -o spectrum_plot.pdf
```

## Clustering mutational spectra

To get an idea of how spectra compare with one another and whether you have groups of similar spectra run (more details in the Clustering_spectra section):

```
python3 cluster_spectra.py -s spectrum1.csv spectrum2.csv spectrum3.csv spectrum4.csv -o output_prefix
```

## Running signal to identify signatures in your spectrum

There's more details on running signal in the Spectrum_decomposition section. The output spectrum from MutTui can be submitted directly to signal (https://signal.mutationalsignatures.com/analyse) as a catalog. If you have multiple spectra you wish to run through signal, you can combine them into a single mutational catalog using:

```
python3 combine_spectra_into_catalog.py -s spectrum1.csv spectrum2.csv -o output_catalog.csv
```

Note that this doesn't combine the spectra into a single spectrum. It creates a file where each of the input spectra is a single column. If you want to combine multiple spectra into a single spectrum, do this first as in the Combining_spectra section.

The spectrum obtained from combine_spectra_into_catalog.py can be submitted to signal as a catalog.
