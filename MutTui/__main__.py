import os, sys
import argparse

from .__init__ import __version__

from .muttui import muttui_parser
from .add_tree_node_labels import add_tree_node_labels_parser
from .cluster_spectra import cluster_spectra_parser
from .combine_spectra import combine_spectra_parser
from .combine_spectra_into_catalog import combine_spectra_into_catalog_parser
from .compare_spectra import compare_spectra_parser
from .convert_vcf_position_translation import convert_vcf_position_translation_parser
from .korimuto import korimuto_parser
from .label_tree import label_tree_parser
from .plot_spectrum import plot_spectrum_parser

def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", title="Available commands")

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    # add subcommands
    muttui_subparser = subparsers.add_parser("run", 
        help="run the MutTui pipeline")
    muttui_subparser = muttui_parser(muttui_subparser)

    korimuto_subparser = subparsers.add_parser("korimuto",
        help="calculates a mutational spectrum from SNP data")
    korimuto_subparser = korimuto_parser(korimuto_subparser)

    plot_spectrum_subparser = subparsers.add_parser("plot",
        help="plot a mutational spectrum")
    plot_spectrum_subparser = plot_spectrum_parser(plot_spectrum_subparser)

    convert_vcf_subparser = subparsers.add_parser("convert-vcf",
        help="creates the position conversion file for MutTui from a VCF file")
    convert_vcf_subparser = convert_vcf_position_translation_parser(convert_vcf_subparser)

    compare_spectra_subparser = subparsers.add_parser("compare",
        help="subtracts one mutational spectrum from another")
    compare_spectra_subparser = compare_spectra_parser(compare_spectra_subparser)

    cluster_spectra_subparser = subparsers.add_parser("cluster", 
        help="cluster spectra using a given tree")
    cluster_spectra_subparser = cluster_spectra_parser(cluster_spectra_subparser)

    combine_spectra_subparser = subparsers.add_parser("combine",
        help="sums the mutations in a number of input spectra to form a single output spectrum")
    combine_spectra_subparser = combine_spectra_parser(combine_spectra_subparser)

    combine_spectra_into_catalog_subparser = subparsers.add_parser("catalogue",
        help="combines multiple input spectra and/or catalogues into a single output catalogue")
    combine_spectra_into_catalog_subparser = combine_spectra_into_catalog_parser(combine_spectra_into_catalog_subparser)

    label_nodes_subparser = subparsers.add_parser("label-nodes", 
        help="add numerical node labels to a given tree to enable labelling of the tree into groups")
    label_nodes_subparser = add_tree_node_labels_parser(label_nodes_subparser)

    label_tree_subparser = subparsers.add_parser("label-tree",
        help="labels nodes in a phylogenetic tree for use with MutTui")
    label_tree_subparser = label_tree_parser(label_tree_subparser)


    # parse arguments and run function
    args = parser.parse_args()
    args.func(args)

    return


if __name__ == "__main__":
    main()
