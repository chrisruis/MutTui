# Calculating spectra from variant data

MutTui can reconstruct mutational spectra from variant (SNP) data as well as from alignments and trees

The input to calculate spectra from variant data can be either a VCF file (default) or a variant file. The variant file is a csv file with columns position (genome position of the variant), reference (containing the reference nucleotide at the position) and variant (containing the variant nucleotide at the position). This file can contain any other columns as long it contains these 3

MutTui identifies the mutations from the VCF or variants file. The contexts are identified from a reference genome. This would normally be the reference genome that was mapped against

To calculate mutational spectra from variant data run:
```
MutTui korimuto -v VCF/variant_file -r reference.fasta -o output_prefix
```

If a variant file is used, option --variant should be specified, for example:
```
MutTui korimuto -v variants.csv -r reference.fasta -o output_prefix --variant
```
