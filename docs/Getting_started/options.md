# What's required

MutTui requires at least 4 files to run:

1) Sequence alignment in fasta format - provided with option -a
2) Phylogenetic tree in newick format - provided with option -t
3) Reference genome in fasta format - provided with option -r
4) Conversion from alignment position to genome position - provided with option -c

The output directory to which the results will be written also need to be specified. The most basic MutTui command is:

```
python3 muttui.py -a alignment.fasta -t tree.nwk -r reference.fasta -c conversion.txt -o muttui_out
```

You can also provide MutTui additional [input files](#dditional-input-files) and [options](#optional-arguments), outlined in the sections below

# Additional input files

#### Labelled tree (-lt)

Use this to give MutTui a labelled tree. MutTui will split each of the labels into a separate mutational spectrum. More details on how to set this up [here](https://github.com/chrisruis/MutTui/tree/main/docs/Labelling_your_tree)

#### Tip labels (-l)

Use this if you want to infer the labels of internal branches from labels at the tips of the tree as part of the MutTui run. MutTui will run the mugration model within treetime to infer internal node labels based on the distribution at tips. This file should contain a header and two columns: the first contains sequence names and the second their label. You don't need to include all sequence names in this file, any sequence names not included are given the label OTHER (or OTHER_A if OTHER is already used in the file)

It's usually best to define your labels before running MutTui where possible and create a labelled tree to provide with -lt, rather than using this -l option

#### Root state (-rs)

Only used if labels are provided with -l and if mugration can't resolve the state at the root. In these rare cases, the root state can be provided

#### GFF (-g)

Use this if you are using either of the options --synonymous or --strand_bias

This should be a GFF file containing gene coordinates in the provided reference genome and is used to identify whether each mutation is synonymous/non-synonymous (if using --synonymous) or which strand each mutation is on (if using --strand_bias)

#### Treetime output directory (-to)

Use this if you have already run treetime ancestral reconstruction and want to run MutTui on those results. This option is only used when --start_from_treetime is specified

MutTui will use the branch_mutations.txt file from the provided directory

#### Minimum branch mutations (-bm)

Use this to set the minimum number of mutations a branch needs to exhibit to have its spectrum calculated. Default is 50

This is only used when --branch_specific is specified which outputs a separate spectrum for each branch with at least -bm mutations

# Optional arguments

#### --rna

**This option is still being developed so doesn't currently work for all scripts**

Use this if working with an RNA pathogen. As it is possible to determine the direction of RNA mutations, reverse mutations are not combined so the output spectrum will contain 192 bars rather than 96

#### --all_sites

Use this if the alignment provided with -a contains all positions in the genome

In this case, the context of each mutation can be identified from the alignment alone so a reference genome does not need to be provided with -r

#### --include_all_branches

By default, when using a labelled tree, the branches on which the label changes are excluded from the calculated spectra as we typically don't know at which point along the branch the label changed. If --include_all_branches is specified, the branches on which the label changes will be included in the spectrum of the downstream label

Use this if specifying labels that only change on tip branches

#### --include_root_branches

Mutations on the two branches immediately downstream of the root are excluded by default as their direction may not be inferred reliably. Use this option to include mutations on these branches

#### --start_from_treetime

Use this option if you have already run treetime ancestral reconstruction and want to use those results for MutTui. The directory containing the treetime output files needs to be provided with -to

#### --strand_bias

Use this option to split the calculated SBS spectrum into transcribed and untranscribed strands, typically to examine strand bias. Mutations are included in the transcribed strand if the pyrimidine is on the transcribed strand and included in the nontranscribed strand if the pyrimidine is on the nontranscribed strand. Mutations that are not within a gene or are in multiple overlapping genes with different orientations are not included in these spectra. Mutations are included if they are in multiple overlapping genes that are in the same orientation

This outputs 2 additional files per branch label - strand_bias_label_X.csv containing the mutation counts and strand_bias_label_X.pdf containing a plot of the mutation counts

#### --synonymous

Use this to only include synonymous mutations in the calculated spectrum. A GFF file needs to be provided with -g

#### --branch_specific

Use this to calculate the output a separate mutational spectrum for each branch in the tree that has at least a given number of mutations. Specify the minimum number of mutations with -bm (default 50)

#### --add_treetime_cmds

Use this to specify additional options that are provided to treetime for ancestral reconstruction

# MutTui usage

```
usage: MutTui run [-h] -a ALIGNMENT -t TREE -o OUTPUT_DIR [-l LABELS] [-lt LABELLED_TREE]
                  [-r REFERENCE] [-c CONVERSION] [-g GFF] [-to TREETIME_OUT] [-bm BRANCH_MUTATIONS]
                  [--add_treetime_cmds ADD_TREETIME_CMDS] [--rna] [--all_sites]
                  [--include_all_branches] [--exclude_root_branches] [-rs ROOT_STATE]
                  [--start_from_treetime] [--strand_bias] [--synonymous] [--branch_specific]
                  [--version]

Run the MutTui pipeline on a given alignment and tree

options:
  -h, --help            show this help message and exit
  --rna                 Specify if using an RNA pathogen, MutTui will output an RNA mutational
                        spectrum
  --all_sites           Specify that the alignment contains all sites, in which case a reference
                        genome does not need to be provided
  --include_all_branches
                        Use when specifying a labelled tree with -lt. By default, the branches along
                        which the label changes are excluded as it is often not clear at what point
                        along the branch the label changed. With adding --include_all_branches, all
                        of the branches along which a label changes will be included in the spectrum
                        of the downstream label
  --exclude_root_branches
                        Use when excluding the two branches that branch directly from the root node.
                        These branches will not be included in any spectrum
  -rs ROOT_STATE, --root_state ROOT_STATE
                        Specify the root state of the given label if known. The root state is
                        calculated by treetime mugration and in almost all cases will be resolved. If
                        the root state has multiple assignments that are equally likely, it cannot be
                        assigned from treetime mugration output. In these rare cases, use this option
                        to assign the root state. This option can also be used to assign a root state
                        if you'd like but its recommended to use the mugration state
  --start_from_treetime
                        Use this option to start with treetime output and so skip inference of
                        ancestral mutations. Use this if you have already run treetime. The directory
                        containing the treetime output files needs to be provided with -to
  --strand_bias         Split mutations into transcribed and untranscribed strands to test for
                        transcription strand bias. A GFF file will need to be provided with option -g
                        to identify genes
  --synonymous          Use non-coding and synonymous mutations only to calculate the mutational
                        spectrum. A GFF file will need to be provided with option -g which will be
                        used to identify genes
  --branch_specific     Calculate the mutational spectrum for each branch in the tree separately.
                        Only branches containing at least the number of mutations specified with -bm
                        will be included
  --version             show program's version number and exit

Input/output:
  -a ALIGNMENT, --alignment ALIGNMENT
                        Input fasta alignment
  -t TREE, --tree TREE  Newick tree file
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        Location of output directory, should already be created and ideally be empty
  -l LABELS, --labels LABELS
                        Comma separated file with sequences and their clade labels to divide the
                        spectrum. If this option is included, treetime mugration will be run to
                        identify changes in the label across the tree. Does not need to include all
                        taxa. Taxa not included will be given a label of OTHER by default, or OTHER_A
                        if OTHER is already taken. This file should have a header.
  -lt LABELLED_TREE, --labelled_tree LABELLED_TREE
                        Tree with state labels, should be the same tree as that provided with -t but
                        with nodes labelled with their state from label_tree.py
  -r REFERENCE, --reference REFERENCE
                        Reference genome sequence containing all sites, used to get context of
                        mutations, not required if using --all_sites
  -c CONVERSION, --conversion CONVERSION
                        Conversion file for alignment position to genome position, used to get
                        context of mutations, not required if using --all_sites
  -g GFF, --gff GFF     GFF reference containing gene coordinates in reference sequence. Used to
                        split mutations into transcription strands and identify synonymous mutations
                        when --synonymous is used
  -to TREETIME_OUT, --treetime_out TREETIME_OUT
                        The location of the directory containing treetime output files from ancestral
                        reconstruction. Only used if option --start_from_treetime is specified, in
                        which case the output files in this directory will be used to calculate the
                        mutational spectrum
  -bm BRANCH_MUTATIONS, --branch_mutations BRANCH_MUTATIONS
                        The minimum number of mutations on a branch to calculate its mutational
                        spectrum. Only used when specifying --branch_specific which will output the
                        spectrum for each branch with more than this number of mutations separately.
                        Default = 50

treetime command:
  --add_treetime_cmds ADD_TREETIME_CMDS
                        Additional options to supply to treetime (these are not checked). Supply
                        these together in quotes
```
