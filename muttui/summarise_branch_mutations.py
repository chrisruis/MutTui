#Calculates the number and proportion of each mutation type on each branch
#Writes a csv summary containing the number and proportion of mutations of each type on each branch
#Writes a labelled tree containing the number and proportion of each mutation on each branch
#Run on all_included_mutations.csv from MutTui and annotated_tree.nexus from treetime (in the MutTui output)
#To run: python3 summarise_branch_mutations.py -m all_included_mutations.csv -t annotated_tree.nexus -o output_prefix

import argparse
from Bio import Phylo
import os

#Extracts depths from a tree, returns a dictionary with clade names as keys and depths as values
def extractDepths(tree):
    #Extract depths of all nodes
    depths = tree.depths()

    #Node names as keys, depths as values
    nodeDepths = dict()

    #Root branch length, often not zero so need to subtract from other depths
    rootBranch = tree.root.branch_length

    #Iterate through the nodes and add to nodeDepths
    for clade in tree.find_clades():
        if clade.is_terminal():
            nodeDepths[clade.name] = depths[clade] - rootBranch
        else:
            nodeDepths[clade.confidence] = depths[clade] - rootBranch
    
    return(nodeDepths)

#Extracts branch lengths from a tree, returns a dictionary with clade names as keys and branch lengths as values
def extractBLs(tree):
    #Node names as keys, branch lengths as values
    bl = dict()

    #Iterate through the nodes and add to bl
    for clade in tree.find_clades():
        if clade.is_terminal():
            bl[clade.name] = clade.branch_length
        else:
            bl[clade.confidence] = clade.branch_length
    
    return(bl)

#Labels each branch in a tree with its total mutations and mutation proportions
def convertTreeMutations(tree, bM, bMT):
    #Iterate through the clades and label them with their mutation proportions
    for clade in tree.find_clades():
        #print(clade.name, depths[clade])
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:

            if clade.name:
                cladeName = clade.name
            else:
                cladeName = clade.confidence
            
            if clade.is_terminal():
                clade.name = clade.name + "[&"
            elif clade.confidence:
                clade.name = "[&label=" + clade.confidence
            else:
                clade.name = "[&label=" + clade.name
            
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
                
                if clade.is_terminal():
                    clade.name += "total_mutations=" + str(totalM)
                else:
                    clade.name += ",total_mutations=" + str(totalM)
                clade.name += ",C_A=" + str(CA) 
                clade.name += ",C_G=" + str(CG)
                clade.name += ",C_T=" + str(CT)
                clade.name += ",T_A=" + str(TA)
                clade.name += ",T_C=" + str(TC)
                clade.name += ",T_G=" + str(TG)
                clade.name += ",P_C_A=" + str(P_CA)
                clade.name += ",P_C_G=" + str(P_CG)
                clade.name += ",P_C_T=" + str(P_CT)
                clade.name += ",P_T_A=" + str(P_TA)
                clade.name += ",P_T_C=" + str(P_TC)
                clade.name += ",P_T_G=" + str(P_TG)
            
            #If the clade has no mutations, all of its mutations will be 0
            else:
                clade.name += ",total_mutations=0,C_A=0,C_G=0,C_T=0,T_A=0,T_C=0,T_G=0,P_C_A=0,P_C_G=0,P_C_T=0,P_T_A=0,P_T_C=0,P_T_G=0"
            
            clade.name += "]"
        
        else:
            clade.name = None
        
        clade.confidence = None
        clade.comment = None
    
    return(tree)

#Convert node labels imported as confidence values to labels
def labelNodeNames(tree):
    for clade in tree.get_nonterminals():
        if clade.confidence:
            clade.label = clade.confidence
        else:
            clade.label = clade.name
    return(tree)

#Converts tip names to a dictionary, used to write trees in BEAST format
def convertTipDict(tree):
    tips = dict()

    tipList = list()

    for i, tip in enumerate(tree.get_terminals()):
        tips[tip.name] = i + 1
        tipList.append(tip.name)
    
    return(tips, tipList)

if __name__ == "__main__":
    description = "Generates a branch-specific summary of mutations in a mutational spectrum"
    parser = argparse.ArgumentParser(description = description)

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

    #Open output files
    outFile = open(args.outfile + ".csv", "w")
    outFile.write("Branch,Branch_depth,Branch_length,Total_mutations,Transitions,Transversions,Tr_Tv_ratio,C>A,C>G,C>T,T>A,T>C,T>G,C>A_proportion,C>G_proportion,C>T_proportion,T>A_proportion,T>C_proportion,T>G_proportion\n")
    outFile_tree = open(args.outfile + "_temp.nex", "w")
    outFile_tree_2 = open(args.outfile + ".nex", "w")

    #Will be filled with branches as keys and total number of mutations as values
    bM = dict()
    #Will be filled with branch:mutation as keys and number of mutations as values
    bMT = dict()

    #Import the tree
    tree = Phylo.read(args.tree.name, "nexus")
    #Extract tree depths
    nodeDepths = extractDepths(tree)
    #Extract tree branch lengths
    branchLengths = extractBLs(tree)

    #The mutation types to be examined
    mutationTypes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    transitions = ["C>T", "T>C"]
    transversions = ["C>A", "C>G", "T>A", "T>G"]

    #Iterate through the mutations and add to bM and bMT
    with open(args.mutations.name) as fileobject:
        next(fileobject)
        for line in fileobject:

            m = line.strip().split(",")[2].split("[")[1].split("]")[0]
            b = line.strip().split(",")[3]

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
        outFile.write(branch + "," + str(nodeDepths[branch]) + "," + str(branchLengths[branch]) + "," + str(bM[branch]) + ",")
        tM = float(bM[branch])

        #Calculate transition/transversion ratio
        tr = float(0)
        tv = float(0)
        for eM in transitions:
            if (branch + ":" + eM) in bMT:
                tr += float(bMT[branch + ":" + eM])
        for eM in transversions:
            if (branch + ":" + eM) in bMT:
                tv += float(bMT[branch + ":" + eM])
        
        outFile.write(str(tr) + "," + str(tv) + ",")
        if tv != float(0):
            outFile.write(str(tr/tv) + ",")
        else:
            outFile.write("NA,")

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

    #Convert the tips so they can be written in BEAST format
    tipDict, tipList = convertTipDict(tree)
    #Change the tip names in the tree to their numeric identifier
    for tip in tree.get_terminals():
        #Add each numeric identifier to bM and bMT so the tip stats can be added to the tree
        if tip.name in bM:
            bM[str(tipDict[tip.name])] = bM[tip.name]
        for mT in mutationTypes:
            if (tip.name + ":" + mT) in bMT:
                bMT[str(tipDict[tip.name]) + ":" + mT] = bMT[tip.name + ":" + mT]
        tip.name = str(tipDict[tip.name])

    #Convert node labels from confidence to labels
    #tree = labelNodeNames(tree)
    #Label each branch with its mutation proportions
    labelledTree = convertTreeMutations(tree, bM, bMT)
    
    #Write the labelled tree
    outFile_tree.write("#NEXUS\n\nBegin Taxa;\n\tDimensions NTax=")
    outFile_tree.write(str(len(tipList)))
    outFile_tree.write(";\n\t\tTaxLabels\n")
    for t in tipList:
        outFile_tree.write("\t\t\t" + t + "\n")
    outFile_tree.write("\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n")
    for t in tipDict:
        outFile_tree.write("\t\t")
        nD = len(str(tipDict[t]))
        if nD == 1:
            outFile_tree.write("   ")
        elif nD == 2:
            outFile_tree.write("  ")
        elif nD == 3:
            outFile_tree.write(" ")
        outFile_tree.write(str(tipDict[t]) + " " + t + ",\n")
    outFile_tree.write(";\ntree TREE1 = ")
    Phylo.write(labelledTree, outFile_tree, "newick")
    outFile_tree.write("End;")
    
    outFile.close()
    outFile_tree.close()

    #Remove quotes from the tree
    treeString = open(outFile_tree.name).readlines()
    for line in treeString:
        outFile_tree_2.write(line.replace("\n", "").replace("'", "") + "\n")
    
    #Remove the temp tree file
    if os.path.exists(args.outfile + "_temp.nex"):
        os.remove(args.outfile + "_temp.nex")

    outFile_tree_2.close()