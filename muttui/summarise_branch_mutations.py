#Calculates the number and proportion of each mutation type on each branch
#Writes a csv summary containing the number and proportion of mutations of each type on each branch
#Writes a labelled tree containing the number and proportion of each mutation on each branch
#Run on all_included_mutations.csv from MutTui and annotated_tree.nexus from treetime (in the MutTui output)
#To run: python3 summarise_branch_mutations.py -m all_included_mutations.csv -t annotated_tree.nexus -o output_prefix

import argparse
from Bio import Phylo
import re

#Takes a string and returns all indices of a given character
def getAllIndices(s, c):
    return([i for i, ltr in enumerate(s) if ltr == c])

#Replaces the mutations in an annotated tree file with mutation proportions
#Takes the tree string from annotated_tree.nexus from treetime
def convertTreeMutations(tree, bM, bMT, mutationTypes):
    #Extract the colon positions in the tree string, used to identify branches
    bPos = getAllIndices(tree, ":")

    #Iterate through the branches, identify their name and replace the mutations with mutation proportions
    for b in bPos[:4]:
        bName = re.split(r"\(|\)|,", tree[:b])[-1]

        if bName in bM:
            totalM = float(bM[bName])
            for m in mutationTypes:
                if (bName + ":" + "C>A") in bMT:
                    CA = float(bMT[bName + ":" + "C>A"])
                    PCA = CA/totalM
                else:
                    CA = "0"
                    PCA = "0"
                if (bName + ":" + "C>G") in bMT:
                    CG = float(bMT[bName + ":" + "C>G"])
                    PCG = CG/totalM
                else:
                    CG = "0"
                    PCG = "0"
                if (bName + ":" + "C>T") in bMT:
                    CT = float(bMT[bName + ":" + "C>T"])
                    PCT = CT/totalM
                else:
                    CT = "0"
                    PCT = "0"
                if (bName + ":" + "T>A") in bMT:
                    TA = float(bMT[bName + ":" + "T>A"])
                    PTA = TA/totalM
                else:
                    TA = "0"
                    PTA = "0"
                if (bName + ":" + "T>C") in bMT:
                    TC = float(bMT[bName + ":" + "T>C"])
                    PTC = TC/totalM
                else:
                    TC = "0"
                    PTC = "0"
                if (bName + ":" + "T>G") in bMT:
                    TG = float(bMT[bName + ":" + "T>G"])
                    PTG = TG/totalM
                else:
                    TG = "0"
                    PTG = "0"
        else:
            totalM = "0"
            CA = "0"
            PCA = "0"
            CG = "0"
            PCG = "0"
            CT = "0"
            PCT = "0"
            TA = "0"
            PTA = "0"
            TC = "0"
            PTC = "0"
            TG = "0"
            PTG = "0"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m",
                        dest = "mutations",
                        required = True,
                        help = "all_included_mutations.csv from MutTui",
                        type = argparse.FileType("r"))
    parser.add_argument("-t",
                        dest = "tree",
                        required = True,
                        help = "annotated_tree.nexus from treetime (in MutTui output directory)",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        dest = "outfile",
                        required = True,
                        help = "Output file prefix. Output csv will be this followed by .csv. " +
                        "Output tree will be this followed by .nex")
    args = parser.parse_args()

    #Open output file
    outFile = open(args.outfile + ".csv", "w")
    outFile.write("Branch,Total_mutations,C>A,C>G,C>T,T>A,T>C,T>G,C>A_proportion,C>G_proportion,C>T_proportion,T>A_proportion,T>C_proportion,T>G_proportion\n")
    outFileTree = open(args.outfile + ".nex", "w")

    #Will be filled with branches as keys and total number of mutations as values
    bM = dict()
    #Will be filled with branch:mutation as keys and number of mutations as values
    bMT = dict()

    #The mutation types to be examined
    mutationTypes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

    #Iterate through the mutations and add to bM and bMT
    with open(args.mutations.name) as fileobject:
        for line in fileobject:

            if line.strip().split(",")[0] != "Mutation_in_alignment":
                m = line.strip().split(",")[2].split("[")[1].split("]")[0]
                b = line.strip().split(",")[3].split("->")[1]

                if b in bM:
                    bM[b] += 1
                else:
                    bM[b] = 1
                
                if (b + ":" + m) in bMT:
                    bMT[b + ":" + m] += 1
                else:
                    bMT[b + ":" + m] = 1
    
    #Write the mutations on each branch
    for branch in bM:
        outFile.write(branch + "," + str(bM[branch]) + ",")
        tM = float(bM[branch])

        for mt in mutationTypes:
            if (branch + ":" + mt) in bMT:
                outFile.write(str(bMT[branch + ":" + mt]) + ",")
            else:
                outFile.write("0,")
        
        for mt in mutationTypes:
            if (branch + ":" + mt) in bMT:
                outFile.write(str(float(bMT[branch + ":" + mt])/tM) + ",")
            else:
                outFile.write("0,")
        
        outFile.write("\n")
    
    #Import the tree as a string
    #tree = Phylo.read(args.tree.name, "nexus")
    tree = open(args.tree.name).readlines()
    #Iterate through the lines in the nexus tree
    #If they are not the tree, write them to the tree outfile
    #If they are the tree, convert the mutations part of the tree to mutation proportions
    for line in tree:
        if line[:5] == " Tree":
            convertTreeMutations(line, bM, bMT, mutationTypes)
        else:
            outFileTree.write(line)
    
    outFile.close()
    outFileTree.close()