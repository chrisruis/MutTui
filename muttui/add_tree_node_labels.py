#Adds node labels to a tree
#First processing step to label a tree for MutTui

import argparse
from Bio import Phylo

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick tree file",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Output tree file")
    args = parser.parse_args()

    #Label the tree nodes
    labelledTree = labelTreeNodes(Phylo.read(args.tree.name, "newick"))

    #Write the labelled tree
    Phylo.write(labelledTree, args.outFile, "newick")