#Labels a tree for use with MutTui based on branches where a given gene mutates

import argparse
import pandas as pd
from Bio import Phylo

#Extract mutations from a variant effect file
#Returns dictionary with node names as keys and mutations as values
def extractMutationEffects(mFile, genes, effects):
    #Import mutations
    mutations = pd.read_csv(mFile.name, sep = "\t")

    #Empty dictionary
    mD = dict()

    #Iterate through the mutations, check if they are in a gene of interest
    #with an effect of interest and if so add to mD
    for i in range(mutations.shape[0]):
        if (mutations["locus_tag"][i] in genes) and (mutations["impact"][i] in effects):
            if mutations["node"][i] in mD:
                mD[mutations["node"][i]].append(mutations["locus_tag"][i])
            else:
                mD[mutations["node"][i]] = [mutations["locus_tag"][i]]
    
    return(mD)

#From label_tree.py
#Gets the state of the upstream clade
def getParentState(tree, clade):
    #Check if the clade is one downstream of the root
    if len(tree.get_path(clade)) >= 2:
        nodeState = tree.get_path(clade)[-2].state
    else:
        nodeState = str(tree.root.state)
    
    return(nodeState)

#Labels a tree based on gene mutations
def labelTreeMutations(tree, mutations):
    #Label the root
    tree.root.state = "R"

    #Keys are branch names, values are labels
    branches = dict()

    #Iterate through the clades in the tree, check if they have mutations in mutations,
    #update their label if necessary, if not add label from parent clade
    for clade in tree.find_clades():
        upstreamLabel = getParentState(tree, clade)
        
        #Add the gene label if there are mutations on the branch
        if clade.name in mutations:
            #Reset the label if it hasn't changed since the root
            if upstreamLabel == "R":
                upstreamLabel = ""
            
            #Add the genes to the node label
            for gene in mutations[clade.name]:
                if gene not in upstreamLabel.split("...."):
                    if upstreamLabel == "":
                        upstreamLabel = gene
                    else:
                        upstreamLabel = upstreamLabel + "...." + gene
            
        clade.state = upstreamLabel
    
    #Relabel the tree names with the state so the tree can be saved with labels
    #Append to the taxon name if a tip branch, assign to the clade name if an internal branch
    for clade in tree.find_clades():
        branches[clade.name] = clade.state
        if clade.is_terminal():
            clade.name = clade.name + "____" + clade.state
        else:
            clade.name = clade.state
    
    return(tree, branches)

if __name__ == "__main__":
    description = "Labels a tree for MutTui based on mutations in one or more genes"
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick tree file to the tree that will be labelled",
                        type = argparse.FileType("r"))
    parser.add_argument("-m",
                        "--mutations",
                        dest = "mutations",
                        required = True,
                        help = "Mutation effect prediction file. This needs to contain columns: node containing branch names that match those " + 
                        "in the tree provided with -t, locus_tag containing the gene names and impact containing the predicted mutation effect",
                        type = argparse.FileType("r"))
    parser.add_argument("-g",
                        "--genes",
                        dest = "genes",
                        required = True,
                        nargs = "+",
                        help = "List of genes that will be used to label the tree. Mutations in any of the " + 
                        "genes in this list with one of the effects provided with -e will be labelled")
    parser.add_argument("-e",
                        "--effect",
                        dest = "effect",
                        nargs = "+",
                        help = "The mutation effects that will result in a mutation being added to the tree. " + 
                        "By default, this will include LOW, MODERATE and HIGH effects so will include all mutations. Provide " + 
                        "one or more of those labels to only include those effects",
                        default = ["LOW", "MODERATE", "HIGH"])
    parser.add_argument("-o",
                        "--outFile",
                        dest = "outFile",
                        help = "Output file prefix. 2 output files are writte: a .nwk file containing the labelled tree " + 
                        "and a .csv file containg the branch in column 1 and its label in column 2")
    
    args = parser.parse_args()

    #Import the tree
    tree = Phylo.read(args.tree.name, "newick")

    #Import the mutations and extract to dictionary with node names as keys and mutations to be labelled as values
    mutations = extractMutationEffects(args.mutations, args.genes, args.effect)

    #Label the tree based on mutations
    stateTree, branches = labelTreeMutations(tree, mutations)

    #Write the tree
    Phylo.write(stateTree, args.outFile + ".nwk", "newick")

    #Write the branch conversion
    outFile = open(args.outFile + ".csv", "w")
    outFile.write("branch,label\n")

    for branch in branches:
        outFile.write(str(branch) + "," + branches[branch] + "\n")

    outFile.close()