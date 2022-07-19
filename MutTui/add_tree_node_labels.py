#Adds node labels to a tree
#First processing step to label a tree for MutTui

import argparse
from Bio import Phylo

#Removes confidence values from the tree
#These will be bootstrap supports, etc that will be written with the node labels if not removed
def cleanTree(tree):
    #Iterate through the clades and remove clade confidence values
    for clade in tree.find_clades():
        if clade.confidence or (clade.confidence == 0.0):
            clade.confidence = None
    
    return(tree)

#Adds node labels to a given tree
def labelTreeNodes(tree):
    #Will be incremented with each node
    nodeIterator = 1

    #Iterate through the internal nodes and add labels
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = "Node" + str(nodeIterator)
            nodeIterator += 1
    
    return(tree)

def add_tree_node_labels_parser(parser):

    parser.description = "Adds numerical node labels to a given tree to enable labelling of the tree into groups"

    parser.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick tree file to be labelled",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Name of output tree file")
    
    parser.set_defaults(func=add_tree_node_labels)
    
    return(parser)


def add_tree_node_labels(args):

    #Clean the tree to remove any bootstrap supports
    tree = cleanTree(Phylo.read(args.tree.name, "newick"))
    
    #Label the tree nodes
    labelledTree = labelTreeNodes(tree)

    #Write the labelled tree
    Phylo.write(labelledTree, args.outFile, "newick")

    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = add_tree_node_labels_parser(parser)
    args = parser.parse_args()

    # run add_tree_node_labels
    args.func(args)

    return

if __name__ == "__main__":
    main()
