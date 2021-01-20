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

#Labels each branch in a tree with its total mutations and mutation proportions
def convertTreeMutations(tree, bM, bMT):
    #Iterate through the clades and label them with their mutation proportions
    for clade in tree.find_clades():
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:

            clade.comment = "&"

            if clade.name:
                cladeName = clade.name
            else:
                cladeName = clade.confidence
            
            #Label the branch with the total number of mutations and mutations of each type
            if cladeName in bM:
                totalM = float(bM[cladeName])

                #Calculate the number and proportion of each mutation type
                if (cladeName + ":C>A") in bMT:
                    CA = float(bMT[cladeName + ":C>A"])
                    P_CA = CA/totalM
                else:
                    CA = float(0)
                    P_CA = float(0)
                
                if (cladeName + ":C>G") in bMT:
                    CG = float(bMT[cladeName + ":C>G"])
                    P_CG = CG/totalM
                else:
                    CG = float(0)
                    P_CG = float(0)
                
                if (cladeName + ":C>T") in bMT:
                    CT = float(bMT[cladeName + ":C>T"])
                    P_CT = CT/totalM
                else:
                    CT = float(0)
                    P_CT = float(0)
                
                if (cladeName + ":T>A") in bMT:
                    TA = float(bMT[cladeName + ":T>A"])
                    P_TA = TA/totalM
                else:
                    TA = float(0)
                    P_TA = float(0)
                
                if (cladeName + ":T>C") in bMT:
                    TC = float(bMT[cladeName + ":T>C"])
                    P_TC = TC/totalM
                else:
                    TC = float(0)
                    P_TC = float(0)
                
                if (cladeName + ":T>G") in bMT:
                    TG = float(bMT[cladeName + ":T>G"])
                    P_TG = TG/totalM
                else:
                    TG = float(0)
                    P_TG = float(0)
                
                clade.comment += "total_mutations=" + str(totalM)
                clade.comment += ",C>A=" + str(CA) 
                clade.comment += ",C>G=" + str(CG)
                clade.comment += ",C>T=" + str(CT)
                clade.comment += ",T>A=" + str(TA)
                clade.comment += ",T>C=" + str(TC)
                clade.comment += ",T>G=" + str(TG)
                clade.comment += ",P_C>A=" + str(P_CA)
                clade.comment += ",P_C>G=" + str(P_CG)
                clade.comment += ",P_C>T=" + str(P_CT)
                clade.comment += ",P_T>A=" + str(P_TA)
                clade.comment += ",P_T>C=" + str(P_TC)
                clade.comment += ",P_T>G=" + str(P_TG)
            
            #If the clade has no mutations, all of its mutations will be 0
            else:
                clade.comment += "total_mutations=0,C>A=0,C>G=0,C>T=0,T>A=0,T>C=0,T>G=0,P_C>A=0,P_C>G=0,P_C>T=0,P_T>A=0,P_T>C=0,P_T>G=0"
        
        clade.confidence = None
    
    return(tree)

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
    tree = Phylo.read(args.tree.name, "nexus")

    #Label each branch with its mutation proportions and write
    labelledTree = convertTreeMutations(tree, bM, bMT)
    Phylo.write(labelledTree, args.outfile + ".nex", "nexus")
    
    outFile.close()