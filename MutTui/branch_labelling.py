#Functions for labelling branches into spectrum categories

import operator
from Bio import Phylo

#Converts taxon labels to a dictionary with labels as keys and taxa as values
def getLabelDict(tree, labels):
    #Import the labels
    taxonLabels = open(labels).readlines()

    #Labels and keys, taxa as values
    labelDict = {}

    #Taxa that have a label, others will be assigned a label of OTHER
    analysedTaxa = []

    #Iterate through the taxa and add to labelDict
    for taxon in taxonLabels[1:]:
        if taxon.strip().split(",")[1] in labelDict:
            labelDict[taxon.strip().split(",")[1]].append(taxon.strip().split(",")[0])
        else:
            labelDict[taxon.strip().split(",")[1]] = [taxon.strip().split(",")[0]]
        analysedTaxa.append(taxon.strip().split(",")[0])
    
    #Check if OTHER has been used as a label, if it has use another OTHER label
    if "OTHER" in labelDict:
        otherLabel = "OTHER_A"
        print("OTHER has been used as a label, using OTHER_A as the label for taxa not included in the label file")
    else:
        otherLabel = "OTHER"
    
    labelDict[otherLabel] = []

    #Iterate through the taxa and add to labelDict if not already included
    for eachTaxon in tree.get_terminals():
        if eachTaxon.name not in analysedTaxa:
            labelDict[otherLabel].append(eachTaxon.name)
    
    return(labelDict)

#Adds the same label to all branches, used when there are no taxon labels provided
def labelAllBranches(tree):
    #Will contain the labels in the tree
    treeLabels = ["A"]

    #Iterate through the branches and add the label A to them
    for clade in tree.find_clades():
        clade.clade_label = "A"
    
    return(tree, treeLabels)

#Labels branches in a given tree with names as in treetime
#Branches are labelled from NODE_0000000 onwards
def labelBranchesTreetime(tree):
    #Appended with each internal node
    iterator = 0

    for clade in tree.find_clades():
        if not clade.is_terminal():
            #Check if the clade already has a name in the input tree, if it does then keep that name
            #If it doesn't, name it with its node number in the tree
            #Therefore if the tree is partially labelled, nodes with labels will keep their label while
            #nodes without labels will be given a node label. This is how treetime deals with partially labelled
            #trees so will enable matching of branch names with the branch_mutations.txt file
            if not clade.name:
                clade.name = "NODE_" + str(iterator).zfill(7)
            iterator += 1
    
    return(tree)

#Writes the label dictionary to a file that can be used by treetime mugration
def writeLabels(labelDict, out_dir):
    labelFile = open(out_dir + "all_taxon_labels.csv", "w")
    labelFile.write("name,group\n")
    for eachLabel in labelDict:
        for sequence in labelDict[eachLabel]:
            labelFile.write(sequence + "," + eachLabel + "\n")
    labelFile.close()

#Identifies the root state from a treetime mugration file
def rootState(tree, confidence, gtr, root_state):
    if root_state:
        rState = root_state
    
    #Identify the root state from confidence
    else:
        #Identify the characters that correspond to each label
        attributes = {}
        for eachLine in gtr:
            if eachLine == "\n":
                break
            elif eachLine.strip() != "Character to attribute mapping:":
                attributes[eachLine.strip().split(": ")[0]] = eachLine.strip().split(": ")[1]
    
        #States and probabilites at the root
        prob = {}
        #Identify the most likely root state from confidence
        for line in confidence:
            if line.strip().split(",")[0] == "NODE_0000000":
                for i in range(1, len(line.strip().split(","))):
                    prob[confidence[0].strip().split(",")[i]] = float(line.strip().split(",")[i])
    
        #State at the root
        rState = attributes[(max(prob.items(), key = operator.itemgetter(1))[0]).strip()]

    #Add the label to the root
    for clade in tree.find_clades():
        #Check if the clade is the root
        if len(tree.get_path(clade)) == 0:
            clade.comment = '[&group="' + rState + '"]'
    
    return(tree)

#Transfers labels from a tree labelled with label_tree.py to the treetime reconstruction tree
def getLabelledTreeLabels(tree, treeFile):
    nodeTree = Phylo.read(treeFile.name, "newick")

    #Will be incremented with each node
    nodeIterator = 1

    #Extract nodeTree to a dictionary with nodes as keys and labels and tips as values
    #Will use this to transfer labels from nodeTree to the main tree by matching the descendent tips
    #Cannot use the same labels as treetime often iterates through the branches in a different order
    nodeTips = {}

    #Iterate through the clades in nodeTree and add to nodeTips
    for clade in nodeTree.find_clades():
        if clade.is_terminal():
            nodeTips["Node" + str(nodeIterator)] = [clade.name.split("____")[1], [clade.name.split("____")[0]]]
        else:
            nodeTips["Node" + str(nodeIterator)] = [clade.name, []]
            #Add the tips in the clade to the empty list
            for tip in clade.get_terminals():
                nodeTips["Node" + str(nodeIterator)][1].append(tip.name.split("____")[0])
        nodeIterator += 1

    #Will be filled with the unique labels from the tree
    treeLabels = []

    #Iterate through the clades, identify the corresponding clade in nodeTree from the descendent tips
    #add the label as clade_label and add unique labels to treeLabels
    for clade in tree.find_clades():
        cladeTips = []
        for tip in clade.get_terminals():
            cladeTips.append(tip.name)
        
        #Identify the corresponding clade in nodeTree
        for eachClade in nodeTips:
            if sorted(nodeTips[eachClade][1]) == sorted(cladeTips):
                clade.clade_label = nodeTips[eachClade][0]
                if nodeTips[eachClade][0] not in treeLabels:
                    treeLabels.append(nodeTips[eachClade][0])
                break
    
    return(tree, treeLabels)

#Label all branches in the tree with their branch name
def labelBranchesNames(tree):
    #Will contain the labels in the tree
    treeLabels = []

    #Iterate through the branches and add the label A to them
    for clade in tree.find_clades():
        clade.clade_label = clade.name
    
    return(tree, treeLabels)

#Iterates through the tree and labels each branch with its category
#Branches with the same category will be placed in the same spectrum, there will be one spectrum per category
def labelBranchesMugration(tree, mugrationTree):
    #Extract the mugration tree to a dictionary with clade names as keys and labels as values
    cladeDict = {}
    for clade in mugrationTree.find_clades():
        if clade.name:
            cladeDict[clade.name] = clade.comment.split('="')[1].split('"')[0]
        else:
            cladeDict[clade.confidence] = clade.comment.split('="')[1].split('"')[0]

    #The labels in the tree
    treeLabels = list(set(cladeDict.values()))

    #Iterate through the clades in the mutation tree and add their label
    for clade in tree.find_clades():
        if clade.name:
            clade.clade_label = cladeDict[clade.name]
        else:
            clade.clade_label = cladeDict[clade.confidence]
    
    return(tree, treeLabels)