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

    #Clean the tree to remove any bootstrap supports
    tree = cleanTree(Phylo.read(args.tree.name, "newick"))
    
    #Label the tree nodes
    labelledTree = labelTreeNodes(tree)

    #Write the labelled tree
    Phylo.write(labelledTree, args.outFile, "newick")