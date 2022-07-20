#Takes the alignment and tree output by treetime and reconstructs the mutational spectrum

import argparse
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from collections import OrderedDict, defaultdict
import operator
import numpy as np
import array
from treetime import *
import re

translation_table = np.array([[[b'K', b'N', b'K', b'N', b'X'],
                               [b'T', b'T', b'T', b'T', b'T'],
                               [b'R', b'S', b'R', b'S', b'X'],
                               [b'I', b'I', b'M', b'I', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'Q', b'H', b'Q', b'H', b'X'],
                               [b'P', b'P', b'P', b'P', b'P'],
                               [b'R', b'R', b'R', b'R', b'R'],
                               [b'L', b'L', b'L', b'L', b'L'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'E', b'D', b'E', b'D', b'X'],
                               [b'A', b'A', b'A', b'A', b'A'],
                               [b'G', b'G', b'G', b'G', b'G'],
                               [b'V', b'V', b'V', b'V', b'V'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'*', b'Y', b'*', b'Y', b'X'],
                               [b'S', b'S', b'S', b'S', b'S'],
                               [b'*', b'C', b'W', b'C', b'X'],
                               [b'L', b'F', b'L', b'F', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']]])

reduce_array = np.full(200, 4)
reduce_array[[65, 97]] = 0
reduce_array[[67, 99]] = 1
reduce_array[[71, 103]] = 2
reduce_array[[84, 116]] = 3

#Converts the positional translation to a dictionary with alignment positions as keys and genome positions as values
def convertTranslation(positionsFile):
    #Import the translation file
    positions = open(positionsFile.name).readlines()

    conversion = {}

    for eachPosition in positions:
        conversion[int(eachPosition.strip().split("\t")[0])] = int(eachPosition.strip().split("\t")[1])
    
    return(conversion)

#Creates a positional translation dictionary with alignment sites as keys and values
#Used when --all_sites is specified so there is no need to convert the sites
def allSitesTranslation(alignment):
    conversion = {}

    for eachPosition in range(len(alignment[0].seq)):
        conversion[int(eachPosition) + 1] = int(eachPosition) + 1
    
    return(conversion)

#Extracts the mutations and positions along a given branch
def getMutations(branch, translation):
    mutations = []
    for mutation in branch.split('="')[1].split('"')[0].split(","):
        mutations.append([mutation[0], int(mutation[1:-1]), translation[int(mutation[1:-1])], mutation[-1]])

    return(mutations)

#Gets the name of the upstream clade
def getParentName(tree, clade):
    #Check if the clade is one downstream of the root
    if len(tree.get_path(clade)) >= 2:
        #Nodes with no mutations along the upstream branch have their IDs imported as name, nodes with mutations as confidence
        if tree.get_path(clade)[-2].confidence:
            node = tree.get_path(clade)[-2].confidence
        else:
            node = tree.get_path(clade)[-2].name
    else:
        node = str(tree.root)
    return(node)

#Identifies the name of the branch
def getBranchName(tree, clade):
    parentNode = getParentName(tree, clade)

    if clade.is_terminal():
        return(parentNode + "->" + clade.name)
    elif clade.confidence:
        return(parentNode + "->" + clade.confidence)
    else:
        return(parentNode + "->" + clade.name)

#Extracts the mutations to a dictionary with branch names as keys and mutations as values
def getBranchDict(tree, positionTranslation):
    branchDict = {}

    for clade in tree.find_clades():
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:
            #Check if there are mutations along the branch
            if (clade.comment) and (clade.comment != '[&mutations=""]'):
                branchDict[getBranchName(tree, clade)] = getMutations(clade.comment, positionTranslation)
            else:
                branchDict[getBranchName(tree, clade)] = "None"
    
    return(branchDict)

#Extracts the mutations in branch_mutations.txt to a dictionary with branch names as keys and mutations as values
def getBranchMutationNexusDict(NexusFile, translation):
    branchDict = defaultdict(list)

    with open(NexusFile, 'r') as infile:
        matches = re.findall("[^,\(\)]+:[.0-9]+\[\&mutations\=\"[,A-Z0-9]*", infile.read())
        for m in matches:
            bname = m.split(':')[0]
            muts = m.split('="')[1].split(",")
            if muts[0] == '': continue
            for mut in muts:
                branchDict[bname].append([mut[0], int(mut[1:-1]), translation[int(mut[1:-1])], mut[-1]])

    return(branchDict)


#Extracts the mutations in branch_mutations.txt to a dictionary with branch names as keys and mutations as values
# def getBranchMutationDict(branchFile, translation):
#     branchDict = {}

#     #Iterate through the mutations and add to branchDict
#     with open(branchFile) as fileobject:
#         #Skip header
#         next(fileobject)

#         for line in fileobject:
#             node = line.strip().split("\t")[0]

#             if node in branchDict:
#                 branchDict[node].append([line.strip().split("\t")[1], int(line.strip().split("\t")[2]), translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]])
#             else:
#                 branchDict[node] = [[line.strip().split("\t")[1], int(line.strip().split("\t")[2]), translation[int(line.strip().split("\t")[2])], line.strip().split("\t")[3]]]
    
#     return(branchDict)


#Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the alignment are assumed
#and the root sequence from the ancestral reconstruction is used
def getReference(reference, all_sites, alignment, positionTranslation):
    if all_sites:
        #Get the reference as the root of the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                referenceSequence = sequence.seq
    
    #Use the provided reference as the reference sequence and substitute in the nucleotides in the inferred root
    #sequence at each variable position. This means the updateReference function uses the root sequence as a
    #starting sequence and updates based on this
    else:
        ref = AlignIO.read(reference.name, "fasta")
        #Extract the sequence of the reference and ensure it is uppercase
        refSeq = ref[0].seq.upper()

        #Reverse the position translation so genome positions are keys and alignment positions are values
        #Can then iterate through the genome positions and check if they are in this
        toChange = {j:i for i, j in positionTranslation.items()}

        #Extract the root sequences from the reconstruction
        for sequence in alignment:
            if sequence.name == "NODE_0000000":
                root_seq = sequence.seq
        
        #Will be the original reference with the root sequence base substituted at each reconstructed position
        referenceSequence = ""

        #Iterate through the length of the reference genome. If the position is in toChange, take the dictionary value position from
        #the root sequence which only contains positions in the translate file. Otherwise, take the position from the reference
        for pos in range(len(refSeq)):
            #Check if the position is variable and needs to be taken from the root
            if (pos + 1) in toChange:
                referenceSequence += root_seq[toChange[pos + 1] - 1]
            else:
                referenceSequence += refSeq[pos]
    
    return(referenceSequence)

#Updates the original reference sequence so it includes the mutations acquired up to the given clade
#Converts the reference to an array, iterates through the upstream branches and updates each of the positions
#that mutated along an upstream branch to the mutated base
#If a position has changed multiple times along the upstream branches, this keeps the most recent change as the
#clades are iterated through from the root through to the most recent upstream branch
####Used when reading in the annotated tree
####def updateReference(tree, clade, branchMutationDict, refSeq):
    ####Convert the reference sequence to an array
    ####referenceArray = array.array("u", refSeq)

    ####Iterate through the upstream branches leading to the node at the start of the current branch
    ####for upstreamClade in tree.get_path(clade)[:-1]:
        ####Check if there are any mutations along the current upstream branch
        ####if branchMutationDict[getBranchName(tree, upstreamClade)] != "None":
            ####Iterate through the previous mutations and update the reference sequence
            ####for eachMutation in branchMutationDict[getBranchName(tree, upstreamClade)]:
                ####referenceArray[eachMutation[2] - 1] = eachMutation[3]
    
    ####return("".join(referenceArray))

#Updates the original reference sequence so it includes the mutations acquired up to the given clade
#Converts the reference to an array, iterates through the upstream branches and updates each of the positions
#that mutated along an upstream branch to the mutated base
#If a position has changed multiple times along the upstream branches, this keeps the most recent change as the
#clades are iterated through from the root through to the most recent upstream branch
def updateReference(tree, clade, branchMutationDict, refSeq):
    #Convert the reference sequence to an array
    referenceArray = array.array("u", refSeq)

    #Iterate through the upstream branches leading to the node at the start of the current branch
    for upstreamClade in tree.get_path(clade)[:-1]:
        #Check if there are any mutations along the current upstream branch
        if upstreamClade.name in branchMutationDict:
            #Iterate through the previous mutations and update the reference sequence
            for eachMutation in branchMutationDict[upstreamClade.name]:
                referenceArray[eachMutation[2] - 1] = eachMutation[3]
    
    return("".join(referenceArray))

#Identifies the label category of a given branch
#By default returns None if the branch is a transition between labels
#If --include_all_branches is used, no branches are labelled None and transition branches are
#labelled with the downstream node's label
def getBranchCategory(tree, clade, include_all_branches):
    #Identify the label of the upstream node
    #Check if the upstream node is the root
    if len(tree.get_path(clade)) >= 2:
        #Check if the label of the upstream clade is the same as the clade
        if clade.clade_label == tree.get_path(clade)[-2].clade_label:
            return(clade.clade_label)
        elif include_all_branches:
            return(clade.clade_label)
        else:
            return(None)
    
    else:
        if clade.clade_label == tree.root.clade_label:
            return(clade.clade_label)
        elif include_all_branches:
            return(clade.clade_label)
        else:
            return(None)

#Removes mutations that are at the start or end of the genome or do not involve 2 nucleotides
#Writes the removed mutations to outMutationsNotUsed
#Splits double substitutions into a separate list
#Returns the filtered mutations as a list
def filterMutations(branchMutations, clade, nucleotides, referenceLength, outMutationsNotUsed):
    positionsToRemove = []
    doubleSubstitutions = []
    #Used to distinguish double substitutions from longer tracts of substitutions
    #positionsToExclude = []
    #doubleSubstitutionPositions = []

    #print("LENGTH:", len(branchMutations))

    #Check if there are double substitutions along the branch and remove these mutations
    for mutation1 in range(len(branchMutations)):
        #Check if the position is the last mutation, if so don't compare it to the next one. Needed
        #so the final mutation on the branch can still be checked for being at the start or end of the genome
        if (mutation1 + 1) != len(branchMutations):
            #Check if the following mutation is at the adjacent genome position
            if branchMutations[mutation1 + 1][2] == (branchMutations[mutation1][2] + 1):
                if mutation1 in doubleSubstitutions:
                    positionsToRemove.append(mutation1 - 1)
                    positionsToRemove.append(mutation1)
                    positionsToRemove.append(mutation1 + 1)
                #positionsToRemove.append(mutation1)
                #positionsToRemove.append(mutation1 + 1)
                #doubleSubstitutions.append(branchMutations[mutation1])
                #doubleSubstitutions.append(branchMutations[mutation1 + 1])
                doubleSubstitutions.append(mutation1)
                doubleSubstitutions.append(mutation1 + 1)
        #Check if the position is at the start or end of the genome
        if (branchMutations[mutation1][2] == 1) or (branchMutations[mutation1][2] == referenceLength):
            positionsToRemove.append(mutation1)
            #Write the mutation to the mutations not used file
            outMutationsNotUsed.write(branchMutations[mutation1][0] + str(branchMutations[mutation1][1]) + branchMutations[mutation1][3] + "," + branchMutations[mutation1][0] + str(branchMutations[mutation1][2]) + branchMutations[mutation1][3] + "," + clade.name + ",End_of_genome\n")
        #Check if the mutation does not involve 2 nucleotides
        elif ((branchMutations[mutation1][0] not in nucleotides) or (branchMutations[mutation1][3] not in nucleotides)) and (mutation1 not in doubleSubstitutions):
            positionsToRemove.append(mutation1)
            outMutationsNotUsed.write(branchMutations[mutation1][0] + str(branchMutations[mutation1][1]) + branchMutations[mutation1][3] + "," + branchMutations[mutation1][0] + str(branchMutations[mutation1][2]) + branchMutations[mutation1][3] + "," + clade.name + ",Mutation_does_not_involve_two_nucleotides\n")
    
    #If there is a tract of 3 or more substitutions at adjacent positions, some of the mutations
    #will be in doubleSubstitutions twice. Extract the unique double substitution positions
    if len(doubleSubstitutions) != 0:
        doubleSubstitutionPositions = [ position for position in doubleSubstitutions if position not in positionsToRemove ]
        doubleSubstitutionMutations = [ branchMutations[i] for i in doubleSubstitutionPositions ]
        #Ensure the double substitutions are sorted by genome position
        doubleSubstitutionMutations.sort(key = lambda x: x[2])
        #Write the positions in a 3 or more substitution tract
        for s in set(doubleSubstitutions):
            if s in positionsToRemove:
                outMutationsNotUsed.write(branchMutations[s][0] + str(branchMutations[s][1]) + branchMutations[s][3] + "," + branchMutations[s][0] + str(branchMutations[s][2]) + branchMutations[s][3] + "," + clade.name + ",In_tract_of_three_or_more_substitutions\n")
    else:
        doubleSubstitutionMutations = []
    
    #uniqueDoubleSubstitutions = []
    #if len(doubleSubstitutions) != 0:
    #if len(doubleSubstitutionMutations) > 0:
    #    dsSet = set()
    #    for ds in doubleSubstitutions:
    #        if tuple(ds) not in dsSet:
    #            uniqueDoubleSubstitutions.append(ds)
    #            dsSet.add(tuple(ds))
    ####Write the double substitutions
    ####Will be removed once double substitutions are incorporated
    #for uds in uniqueDoubleSubstitutions:
    #    outMutationsNotUsed.write(uds[0] + str(uds[1]) + uds[3] + "," + uds[0] + str(uds[2]) + uds[3] + "," + clade.name + ",Double_substitution\n")
    
    #Add the double substitutions into the positions to be removed from the single substitutions
    for eachPosition in doubleSubstitutions:
        if eachPosition not in positionsToRemove:
            positionsToRemove.append(eachPosition)
    positionsToRemove = set(positionsToRemove)

    #Remove the positions that will not be included in the spectrum
    if len(positionsToRemove) != 0:
        #Identify the unique set of positions, if there are 3 or more consecutive positions to be removed,
        #at least one of these positions will be in positionsToRemove more than once
        uniquePositionsToRemove = list(set(positionsToRemove))
        #Remove the positions from the mutations
        for ele in sorted(uniquePositionsToRemove, reverse = True):
            del branchMutations[ele]
    
    return(branchMutations, doubleSubstitutionMutations)

#Translates a given nucleotide sequence to protein
def translateSequence(sequence, strand):
    if strand == "+":
        return(Seq("".join(sequence)).translate())
    else:
        return(Seq("".join(sequence)).reverse_complement().translate())

#Identifies the amino acid position of a nucleotide within a gene
#If the gene is on the positive strand, this is the nucleotide position divided by 3, rounded up
#If the gene is on the negative strand, this is the nucleotide position from the end of the gene divided by 3, rounded up
def extractPosition(geneCoordinates, positionInGene):
    if geneCoordinates[2] == "+":
        return(int(positionInGene/3) + (positionInGene % 3 > 0))
    else:
        reversePosition = geneCoordinates[1] - (geneCoordinates[0] + positionInGene) + 1
        return(int(reversePosition/3) + (reversePosition % 3 > 0))

#Extracts synonymous mutations along a given branch
#Used when --synonymous is specified
def extractSynonymous(branchMutations, updatedReference, geneCoordinates, positionGene):
    #Gene sequences at the upstream node
    upstreamGenes = dict()
    #Gene sequences containing mutations along the branch
    downstreamGenes = dict()

    positionsToRemove = []

    #Iterate through the mutations to create dictionaries of each gene that mutates along the branch
    #with their nucleotide sequences pre and post mutations
    for mutation in branchMutations:
        #Check if the mutation is in a gene, if the position isn't in positionGene, the position is intergenic
        if mutation[2] not in positionGene:
            continue
        #Iterate through the genes the position is in, check if they are already present, mutate if so, if not
        #add to dictionaries and mutate
        for geneName in positionGene[mutation[2]].split("____"):
            #Position of the mutation in the gene, zero based
            positionInGene = mutation[2] - geneCoordinates[geneName][0]
            if geneName in downstreamGenes:
                downstreamGenes[geneName][positionInGene] = mutation[3]
            else:
                upstreamGenes[geneName] = array.array("u", updatedReference[(geneCoordinates[geneName][0] - 1):geneCoordinates[geneName][1]])
                downstreamGenes[geneName] = array.array("u", updatedReference[(geneCoordinates[geneName][0] - 1):geneCoordinates[geneName][1]])
                #Mutate the position in the downstream gene
                downstreamGenes[geneName][positionInGene] = mutation[3]
    
    #Iterate through the mutations again, check if the translated position of the mutation changes along the branch
    #Second iteration needed because 1st and 3rd codon positions may both change along a branch
    for i, mutation in enumerate(branchMutations):
        #Check if the mutation is in a gene, if the position isn't in positionGene, the position is intergenic
        if mutation[2] not in positionGene:
            continue
        #Will change to False if the amino acid of the mutation has changed in any gene it is present in
        synonymous = True
        #Iterate through the genes the position is in, check if the position's amino acid changes, if so change synonymous to False
        for geneName in positionGene[mutation[2]].split("____"):
            #Position of the mutation in the gene, zero based
            positionInGene = mutation[2] - geneCoordinates[geneName][0]
            upstreamAA = translateSequence(upstreamGenes[geneName], geneCoordinates[geneName][2])
            downstreamAA = translateSequence(downstreamGenes[geneName], geneCoordinates[geneName][2])
            aaPosition = extractPosition(geneCoordinates[geneName], positionInGene) - 1

            if upstreamAA[aaPosition] != downstreamAA[aaPosition]:
                synonymous = False
        
        #If the mutation is nonsynonymous within any gene, remove it
        if synonymous == False:
            positionsToRemove.append(i)
    
    #Remove the positions that will not be included in the spectrum
    if len(positionsToRemove) != 0:
        #Identify the unique set of positions
        uniquePositionsToRemove = list(set(positionsToRemove))
        #Remove the positions from the mutations
        for ele in sorted(uniquePositionsToRemove, reverse = True):
            del branchMutations[ele]
    
    return(branchMutations)

#Identifies the strand bias of each mutation
def getStrandBias(mutation, updatedReference, geneCoordinates, positionGene):
    #print(mutation)
    #Strands of the genes the mutation is in, used to exclude mutations in multiple genes on opposing strands
    ss = set()

    #Get the strand of each gene the mutation is in and add to ss
    for g in positionGene[mutation[2]].split("____"):
        ss.add(geneCoordinates[g][2])
    
    #Transition mutations
    transitions = ["C", "T"]
    
    #Check if the mutation is in multiple genes on different strands, if so exclude
    if len(ss) == 1:
        if (("+" in ss) and (mutation[0] in transitions)) or (("-" in ss) and (mutation[0] not in transitions)):
            return("t")
        else:
            return("u")
    else:
        return(None)

#Identifies the context of a mutation
def getContext(mutation, updatedReference):
    return(updatedReference[mutation[2] - 2], updatedReference[mutation[2]])

#Takes a base and returns the complement base
def complement(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"

#Creates an empty mutational spectrum dictionary for DNA datasets, i.e. combining symmetric mutations
def getMutationDict():
    mutation = OrderedDict()
    mutation["ACAA"] = 0
    mutation["ACAC"] = 0
    mutation["ACAG"] = 0
    mutation["ACAT"] = 0
    mutation["CCAA"] = 0
    mutation["CCAC"] = 0
    mutation["CCAG"] = 0
    mutation["CCAT"] = 0
    mutation["GCAA"] = 0
    mutation["GCAC"] = 0
    mutation["GCAG"] = 0
    mutation["GCAT"] = 0
    mutation["TCAA"] = 0
    mutation["TCAC"] = 0
    mutation["TCAG"] = 0
    mutation["TCAT"] = 0
    mutation["ACGA"] = 0
    mutation["ACGC"] = 0
    mutation["ACGG"] = 0
    mutation["ACGT"] = 0
    mutation["CCGA"] = 0
    mutation["CCGC"] = 0
    mutation["CCGG"] = 0
    mutation["CCGT"] = 0
    mutation["GCGA"] = 0
    mutation["GCGC"] = 0
    mutation["GCGG"] = 0
    mutation["GCGT"] = 0
    mutation["TCGA"] = 0
    mutation["TCGC"] = 0
    mutation["TCGG"] = 0
    mutation["TCGT"] = 0
    mutation["ACTA"] = 0
    mutation["ACTC"] = 0
    mutation["ACTG"] = 0
    mutation["ACTT"] = 0
    mutation["CCTA"] = 0
    mutation["CCTC"] = 0
    mutation["CCTG"] = 0
    mutation["CCTT"] = 0
    mutation["GCTA"] = 0
    mutation["GCTC"] = 0
    mutation["GCTG"] = 0
    mutation["GCTT"] = 0
    mutation["TCTA"] = 0
    mutation["TCTC"] = 0
    mutation["TCTG"] = 0
    mutation["TCTT"] = 0
    mutation["ATAA"] = 0
    mutation["ATAC"] = 0
    mutation["ATAG"] = 0
    mutation["ATAT"] = 0
    mutation["CTAA"] = 0
    mutation["CTAC"] = 0
    mutation["CTAG"] = 0
    mutation["CTAT"] = 0
    mutation["GTAA"] = 0
    mutation["GTAC"] = 0
    mutation["GTAG"] = 0
    mutation["GTAT"] = 0
    mutation["TTAA"] = 0
    mutation["TTAC"] = 0
    mutation["TTAG"] = 0
    mutation["TTAT"] = 0
    mutation["ATCA"] = 0
    mutation["ATCC"] = 0
    mutation["ATCG"] = 0
    mutation["ATCT"] = 0
    mutation["CTCA"] = 0
    mutation["CTCC"] = 0
    mutation["CTCG"] = 0
    mutation["CTCT"] = 0
    mutation["GTCA"] = 0
    mutation["GTCC"] = 0
    mutation["GTCG"] = 0
    mutation["GTCT"] = 0
    mutation["TTCA"] = 0
    mutation["TTCC"] = 0
    mutation["TTCG"] = 0
    mutation["TTCT"] = 0
    mutation["ATGA"] = 0
    mutation["ATGC"] = 0
    mutation["ATGG"] = 0
    mutation["ATGT"] = 0
    mutation["CTGA"] = 0
    mutation["CTGC"] = 0
    mutation["CTGG"] = 0
    mutation["CTGT"] = 0
    mutation["GTGA"] = 0
    mutation["GTGC"] = 0
    mutation["GTGG"] = 0
    mutation["GTGT"] = 0
    mutation["TTGA"] = 0
    mutation["TTGC"] = 0
    mutation["TTGG"] = 0
    mutation["TTGT"] = 0

    return(mutation)

#Creates an empty strand bias spectrum dictionary for DNA datasets, each mutation is in the dictionary twice:
#once on the transcribed strand (preceded by t) and once on the untranscribed strand (preceded by u)
def getStrandBiasDict():
    mutation = OrderedDict()
    mutation["tACAA"] = 0
    mutation["tACAC"] = 0
    mutation["tACAG"] = 0
    mutation["tACAT"] = 0
    mutation["tCCAA"] = 0
    mutation["tCCAC"] = 0
    mutation["tCCAG"] = 0
    mutation["tCCAT"] = 0
    mutation["tGCAA"] = 0
    mutation["tGCAC"] = 0
    mutation["tGCAG"] = 0
    mutation["tGCAT"] = 0
    mutation["tTCAA"] = 0
    mutation["tTCAC"] = 0
    mutation["tTCAG"] = 0
    mutation["tTCAT"] = 0
    mutation["tACGA"] = 0
    mutation["tACGC"] = 0
    mutation["tACGG"] = 0
    mutation["tACGT"] = 0
    mutation["tCCGA"] = 0
    mutation["tCCGC"] = 0
    mutation["tCCGG"] = 0
    mutation["tCCGT"] = 0
    mutation["tGCGA"] = 0
    mutation["tGCGC"] = 0
    mutation["tGCGG"] = 0
    mutation["tGCGT"] = 0
    mutation["tTCGA"] = 0
    mutation["tTCGC"] = 0
    mutation["tTCGG"] = 0
    mutation["tTCGT"] = 0
    mutation["tACTA"] = 0
    mutation["tACTC"] = 0
    mutation["tACTG"] = 0
    mutation["tACTT"] = 0
    mutation["tCCTA"] = 0
    mutation["tCCTC"] = 0
    mutation["tCCTG"] = 0
    mutation["tCCTT"] = 0
    mutation["tGCTA"] = 0
    mutation["tGCTC"] = 0
    mutation["tGCTG"] = 0
    mutation["tGCTT"] = 0
    mutation["tTCTA"] = 0
    mutation["tTCTC"] = 0
    mutation["tTCTG"] = 0
    mutation["tTCTT"] = 0
    mutation["tATAA"] = 0
    mutation["tATAC"] = 0
    mutation["tATAG"] = 0
    mutation["tATAT"] = 0
    mutation["tCTAA"] = 0
    mutation["tCTAC"] = 0
    mutation["tCTAG"] = 0
    mutation["tCTAT"] = 0
    mutation["tGTAA"] = 0
    mutation["tGTAC"] = 0
    mutation["tGTAG"] = 0
    mutation["tGTAT"] = 0
    mutation["tTTAA"] = 0
    mutation["tTTAC"] = 0
    mutation["tTTAG"] = 0
    mutation["tTTAT"] = 0
    mutation["tATCA"] = 0
    mutation["tATCC"] = 0
    mutation["tATCG"] = 0
    mutation["tATCT"] = 0
    mutation["tCTCA"] = 0
    mutation["tCTCC"] = 0
    mutation["tCTCG"] = 0
    mutation["tCTCT"] = 0
    mutation["tGTCA"] = 0
    mutation["tGTCC"] = 0
    mutation["tGTCG"] = 0
    mutation["tGTCT"] = 0
    mutation["tTTCA"] = 0
    mutation["tTTCC"] = 0
    mutation["tTTCG"] = 0
    mutation["tTTCT"] = 0
    mutation["tATGA"] = 0
    mutation["tATGC"] = 0
    mutation["tATGG"] = 0
    mutation["tATGT"] = 0
    mutation["tCTGA"] = 0
    mutation["tCTGC"] = 0
    mutation["tCTGG"] = 0
    mutation["tCTGT"] = 0
    mutation["tGTGA"] = 0
    mutation["tGTGC"] = 0
    mutation["tGTGG"] = 0
    mutation["tGTGT"] = 0
    mutation["tTTGA"] = 0
    mutation["tTTGC"] = 0
    mutation["tTTGG"] = 0
    mutation["tTTGT"] = 0
    mutation["uACAA"] = 0
    mutation["uACAC"] = 0
    mutation["uACAG"] = 0
    mutation["uACAT"] = 0
    mutation["uCCAA"] = 0
    mutation["uCCAC"] = 0
    mutation["uCCAG"] = 0
    mutation["uCCAT"] = 0
    mutation["uGCAA"] = 0
    mutation["uGCAC"] = 0
    mutation["uGCAG"] = 0
    mutation["uGCAT"] = 0
    mutation["uTCAA"] = 0
    mutation["uTCAC"] = 0
    mutation["uTCAG"] = 0
    mutation["uTCAT"] = 0
    mutation["uACGA"] = 0
    mutation["uACGC"] = 0
    mutation["uACGG"] = 0
    mutation["uACGT"] = 0
    mutation["uCCGA"] = 0
    mutation["uCCGC"] = 0
    mutation["uCCGG"] = 0
    mutation["uCCGT"] = 0
    mutation["uGCGA"] = 0
    mutation["uGCGC"] = 0
    mutation["uGCGG"] = 0
    mutation["uGCGT"] = 0
    mutation["uTCGA"] = 0
    mutation["uTCGC"] = 0
    mutation["uTCGG"] = 0
    mutation["uTCGT"] = 0
    mutation["uACTA"] = 0
    mutation["uACTC"] = 0
    mutation["uACTG"] = 0
    mutation["uACTT"] = 0
    mutation["uCCTA"] = 0
    mutation["uCCTC"] = 0
    mutation["uCCTG"] = 0
    mutation["uCCTT"] = 0
    mutation["uGCTA"] = 0
    mutation["uGCTC"] = 0
    mutation["uGCTG"] = 0
    mutation["uGCTT"] = 0
    mutation["uTCTA"] = 0
    mutation["uTCTC"] = 0
    mutation["uTCTG"] = 0
    mutation["uTCTT"] = 0
    mutation["uATAA"] = 0
    mutation["uATAC"] = 0
    mutation["uATAG"] = 0
    mutation["uATAT"] = 0
    mutation["uCTAA"] = 0
    mutation["uCTAC"] = 0
    mutation["uCTAG"] = 0
    mutation["uCTAT"] = 0
    mutation["uGTAA"] = 0
    mutation["uGTAC"] = 0
    mutation["uGTAG"] = 0
    mutation["uGTAT"] = 0
    mutation["uTTAA"] = 0
    mutation["uTTAC"] = 0
    mutation["uTTAG"] = 0
    mutation["uTTAT"] = 0
    mutation["uATCA"] = 0
    mutation["uATCC"] = 0
    mutation["uATCG"] = 0
    mutation["uATCT"] = 0
    mutation["uCTCA"] = 0
    mutation["uCTCC"] = 0
    mutation["uCTCG"] = 0
    mutation["uCTCT"] = 0
    mutation["uGTCA"] = 0
    mutation["uGTCC"] = 0
    mutation["uGTCG"] = 0
    mutation["uGTCT"] = 0
    mutation["uTTCA"] = 0
    mutation["uTTCC"] = 0
    mutation["uTTCG"] = 0
    mutation["uTTCT"] = 0
    mutation["uATGA"] = 0
    mutation["uATGC"] = 0
    mutation["uATGG"] = 0
    mutation["uATGT"] = 0
    mutation["uCTGA"] = 0
    mutation["uCTGC"] = 0
    mutation["uCTGG"] = 0
    mutation["uCTGT"] = 0
    mutation["uGTGA"] = 0
    mutation["uGTGC"] = 0
    mutation["uGTGG"] = 0
    mutation["uGTGT"] = 0
    mutation["uTTGA"] = 0
    mutation["uTTGC"] = 0
    mutation["uTTGG"] = 0
    mutation["uTTGT"] = 0

    return(mutation)

#Creates an empty mutational spectrum dictionary for RNA datasets, i.e. not combining symmetric mutations
def getRNADict():
    mutation = OrderedDict()
    mutation["ACAA"] = 0
    mutation["ACAC"] = 0
    mutation["ACAG"] = 0
    mutation["ACAT"] = 0
    mutation["CCAA"] = 0
    mutation["CCAC"] = 0
    mutation["CCAG"] = 0
    mutation["CCAT"] = 0
    mutation["GCAA"] = 0
    mutation["GCAC"] = 0
    mutation["GCAG"] = 0
    mutation["GCAT"] = 0
    mutation["TCAA"] = 0
    mutation["TCAC"] = 0
    mutation["TCAG"] = 0
    mutation["TCAT"] = 0
    mutation["ACGA"] = 0
    mutation["ACGC"] = 0
    mutation["ACGG"] = 0
    mutation["ACGT"] = 0
    mutation["CCGA"] = 0
    mutation["CCGC"] = 0
    mutation["CCGG"] = 0
    mutation["CCGT"] = 0
    mutation["GCGA"] = 0
    mutation["GCGC"] = 0
    mutation["GCGG"] = 0
    mutation["GCGT"] = 0
    mutation["TCGA"] = 0
    mutation["TCGC"] = 0
    mutation["TCGG"] = 0
    mutation["TCGT"] = 0
    mutation["ACTA"] = 0
    mutation["ACTC"] = 0
    mutation["ACTG"] = 0
    mutation["ACTT"] = 0
    mutation["CCTA"] = 0
    mutation["CCTC"] = 0
    mutation["CCTG"] = 0
    mutation["CCTT"] = 0
    mutation["GCTA"] = 0
    mutation["GCTC"] = 0
    mutation["GCTG"] = 0
    mutation["GCTT"] = 0
    mutation["TCTA"] = 0
    mutation["TCTC"] = 0
    mutation["TCTG"] = 0
    mutation["TCTT"] = 0
    mutation["ATAA"] = 0
    mutation["ATAC"] = 0
    mutation["ATAG"] = 0
    mutation["ATAT"] = 0
    mutation["CTAA"] = 0
    mutation["CTAC"] = 0
    mutation["CTAG"] = 0
    mutation["CTAT"] = 0
    mutation["GTAA"] = 0
    mutation["GTAC"] = 0
    mutation["GTAG"] = 0
    mutation["GTAT"] = 0
    mutation["TTAA"] = 0
    mutation["TTAC"] = 0
    mutation["TTAG"] = 0
    mutation["TTAT"] = 0
    mutation["ATCA"] = 0
    mutation["ATCC"] = 0
    mutation["ATCG"] = 0
    mutation["ATCT"] = 0
    mutation["CTCA"] = 0
    mutation["CTCC"] = 0
    mutation["CTCG"] = 0
    mutation["CTCT"] = 0
    mutation["GTCA"] = 0
    mutation["GTCC"] = 0
    mutation["GTCG"] = 0
    mutation["GTCT"] = 0
    mutation["TTCA"] = 0
    mutation["TTCC"] = 0
    mutation["TTCG"] = 0
    mutation["TTCT"] = 0
    mutation["ATGA"] = 0
    mutation["ATGC"] = 0
    mutation["ATGG"] = 0
    mutation["ATGT"] = 0
    mutation["CTGA"] = 0
    mutation["CTGC"] = 0
    mutation["CTGG"] = 0
    mutation["CTGT"] = 0
    mutation["GTGA"] = 0
    mutation["GTGC"] = 0
    mutation["GTGG"] = 0
    mutation["GTGT"] = 0
    mutation["TTGA"] = 0
    mutation["TTGC"] = 0
    mutation["TTGG"] = 0
    mutation["TTGT"] = 0
    mutation["AGTA"] = 0
    mutation["AGTC"] = 0
    mutation["AGTG"] = 0
    mutation["AGTT"] = 0
    mutation["CGTA"] = 0
    mutation["CGTC"] = 0
    mutation["CGTG"] = 0
    mutation["CGTT"] = 0
    mutation["GGTA"] = 0
    mutation["GGTC"] = 0
    mutation["GGTG"] = 0
    mutation["GGTT"] = 0
    mutation["TGTA"] = 0
    mutation["TGTC"] = 0
    mutation["TGTG"] = 0
    mutation["TGTT"] = 0
    mutation["AGCA"] = 0
    mutation["AGCC"] = 0
    mutation["AGCG"] = 0
    mutation["AGCT"] = 0
    mutation["CGCA"] = 0
    mutation["CGCC"] = 0
    mutation["CGCG"] = 0
    mutation["CGCT"] = 0
    mutation["GGCA"] = 0
    mutation["GGCC"] = 0
    mutation["GGCG"] = 0
    mutation["GGCT"] = 0
    mutation["TGCA"] = 0
    mutation["TGCC"] = 0
    mutation["TGCG"] = 0
    mutation["TGCT"] = 0
    mutation["AGAA"] = 0
    mutation["AGAC"] = 0
    mutation["AGAG"] = 0
    mutation["AGAT"] = 0
    mutation["CGAA"] = 0
    mutation["CGAC"] = 0
    mutation["CGAG"] = 0
    mutation["CGAT"] = 0
    mutation["GGAA"] = 0
    mutation["GGAC"] = 0
    mutation["GGAG"] = 0
    mutation["GGAT"] = 0
    mutation["TGAA"] = 0
    mutation["TGAC"] = 0
    mutation["TGAG"] = 0
    mutation["TGAT"] = 0
    mutation["AATA"] = 0
    mutation["AATC"] = 0
    mutation["AATG"] = 0
    mutation["AATT"] = 0
    mutation["CATA"] = 0
    mutation["CATC"] = 0
    mutation["CATG"] = 0
    mutation["CATT"] = 0
    mutation["GATA"] = 0
    mutation["GATC"] = 0
    mutation["GATG"] = 0
    mutation["GATT"] = 0
    mutation["TATA"] = 0
    mutation["TATC"] = 0
    mutation["TATG"] = 0
    mutation["TATT"] = 0
    mutation["AAGA"] = 0
    mutation["AAGC"] = 0
    mutation["AAGG"] = 0
    mutation["AAGT"] = 0
    mutation["CAGA"] = 0
    mutation["CAGC"] = 0
    mutation["CAGG"] = 0
    mutation["CAGT"] = 0
    mutation["GAGA"] = 0
    mutation["GAGC"] = 0
    mutation["GAGG"] = 0
    mutation["GAGT"] = 0
    mutation["TAGA"] = 0
    mutation["TAGC"] = 0
    mutation["TAGG"] = 0
    mutation["TAGT"] = 0
    mutation["AACA"] = 0
    mutation["AACC"] = 0
    mutation["AACG"] = 0
    mutation["AACT"] = 0
    mutation["CACA"] = 0
    mutation["CACC"] = 0
    mutation["CACG"] = 0
    mutation["CACT"] = 0
    mutation["GACA"] = 0
    mutation["GACC"] = 0
    mutation["GACG"] = 0
    mutation["GACT"] = 0
    mutation["TACA"] = 0
    mutation["TACC"] = 0
    mutation["TACG"] = 0
    mutation["TACT"] = 0

    return(mutation)

#Creates an empty double substitution mutational spcetrum dictionary for DNA datasets, i.e. combining symmetric mutations
#Each entry is reference 1, reference 2, mutation 1, mutation 2
def getDoubleSubstitutionDict():
    mutation = OrderedDict()
    mutation["AACC"] = 0
    mutation["AACG"] = 0
    mutation["AACT"] = 0
    mutation["AAGC"] = 0
    mutation["AAGG"] = 0
    mutation["AAGT"] = 0
    mutation["AATC"] = 0
    mutation["AATG"] = 0
    mutation["AATT"] = 0
    mutation["ACCA"] = 0
    mutation["ACCG"] = 0
    mutation["ACCT"] = 0
    mutation["ACGA"] = 0
    mutation["ACGG"] = 0
    mutation["ACGT"] = 0
    mutation["ACTA"] = 0
    mutation["ACTG"] = 0
    mutation["ACTT"] = 0
    mutation["AGCA"] = 0
    mutation["AGCC"] = 0
    mutation["AGCT"] = 0
    mutation["AGGA"] = 0
    mutation["AGGC"] = 0
    mutation["AGGT"] = 0
    mutation["AGTA"] = 0
    mutation["AGTC"] = 0
    mutation["AGTT"] = 0
    mutation["ATCA"] = 0
    mutation["ATCC"] = 0
    mutation["ATCG"] = 0
    mutation["ATGA"] = 0
    mutation["ATGC"] = 0
    mutation["ATTA"] = 0
    mutation["CAAC"] = 0
    mutation["CAAG"] = 0
    mutation["CAAT"] = 0
    mutation["CAGC"] = 0
    mutation["CAGG"] = 0
    mutation["CAGT"] = 0
    mutation["CATC"] = 0
    mutation["CATG"] = 0
    mutation["CATT"] = 0
    mutation["CCAA"] = 0
    mutation["CCAG"] = 0
    mutation["CCAT"] = 0
    mutation["CCGA"] = 0
    mutation["CCGG"] = 0
    mutation["CCGT"] = 0
    mutation["CCTA"] = 0
    mutation["CCTG"] = 0
    mutation["CCTT"] = 0
    mutation["CGAA"] = 0
    mutation["CGAC"] = 0
    mutation["CGAT"] = 0
    mutation["CGGA"] = 0
    mutation["CGGC"] = 0
    mutation["CGTA"] = 0
    mutation["GAAC"] = 0
    mutation["GAAG"] = 0
    mutation["GAAT"] = 0
    mutation["GACC"] = 0
    mutation["GACG"] = 0
    mutation["GACT"] = 0
    mutation["GATC"] = 0
    mutation["GATG"] = 0
    mutation["GATT"] = 0
    mutation["GCAA"] = 0
    mutation["GCAG"] = 0
    mutation["GCAT"] = 0
    mutation["GCCA"] = 0
    mutation["GCCG"] = 0
    mutation["GCTA"] = 0
    mutation["TAAC"] = 0
    mutation["TAAG"] = 0
    mutation["TAAT"] = 0
    mutation["TACC"] = 0
    mutation["TACG"] = 0
    mutation["TAGC"] = 0
    return(mutation)

#Calculates the number of each mutation type in the spectrum
#Takes a mutational spectrum dictionary
#Returns a dictionary with mutation types as keys and counts as values
def mutationTypeCount(spectrum, rna):
    #Empty dictionary of mutation counts
    if rna:
        mtDict = {"AC": 0, "AG": 0, "AT": 0, "CA": 0, "CG": 0, "CT": 0, "GA": 0, "GC": 0, "GT": 0, "TA": 0, "TC": 0, "TG": 0}
    else:
        mtDict = {"CA": 0, "CG": 0, "CT": 0, "TA": 0, "TC": 0, "TG": 0}
    
    #Iterate through the mutations in the spectrum, extract their mutation type and add to the corresponding dictionary mutation
    for m in spectrum:
        mtDict[m[1] + m[2]] += spectrum[m]
    
    return(mtDict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "ancestral_sequences.fasta from treetime")
    parser.add_argument("-t", help = "annotated_tree.nexus from treetime")
    parser.add_argument("-r", help = "Reference genome in fasta format")
    parser.add_argument("-p", help = "Translation of alignment positions to genome positions, tab separated, no header, alignment positions in column 1, genome positions in column 2")
    parser.add_argument("-o", help = "Output file prefix, default reconstructions if none given", default = "reconstructions")
    args = parser.parse_args()

    #Import the ancestral sequence alignment
    alignment = AlignIO.read(args.a, "fasta")

    #Import the annotated tree
    tree = Phylo.read(args.t, "nexus")

    #Import the reference genome sequence
    reference = AlignIO.read(args.r, "fasta")
    referenceLength = len(reference[0])
    referenceSequence = str(reference[0].seq)

    #Import the translation file
    alignmentGenome = open(args.p).readlines()

    #Will contain the number of each contextual mutation
    outMutations = open(args.o + "_mutational_spectrum.txt", "w")
    outMutations.write("Mutation_with_context,Number_of_mutations\n")

    #Will contain the mutations that were not included in the contextual spectrum due to lack of surrounding context, includes mutations next to gaps and mutations at the end of the sequence
    outMutationsNotUsed = open(args.o + "_mutations_not_included.txt", "w")
    outMutationsNotUsed.write("Mutation,Branch,Reason_not_included\n")
    
    #Convert the translation file to a dictionary
    positionTranslation = convertTranslation(alignmentGenome)

    #Identify the variable positions in the genome
    variablePositions = getVariablePositions(positionTranslation)

    #Extract the mutations in each branch into a dictionary
    branchMutationDict = getBranchDict(tree, positionTranslation)

    #Create an empty mutational spectrum
    spectrum = getMutationDict()

    #The 4 nucleotides, used to check if upstream and downstream bases are nucleotides
    nucleotides = ["A","C","G","T"]

    #Iterate through the branches, extract the mutations, get their context and add to mutationSpectrum
    for clade in tree.find_clades():
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:
            #Identify the name of the branch for saving mutations
            branchName = getBranchName(tree, clade)
            #Check if there are mutations along the branch, the .comment is only added to the clade if there are
            if clade.comment:
                branchMutations = getMutations(clade.comment, positionTranslation)

                #Check if there are double substitutions along the branch and remove these mutations
                if len(branchMutations) > 1:
                    positionsToRemove = []
                    for mutation1 in range(0, (len(branchMutations) - 1)):
                        for mutation2 in range((mutation1 + 1), len(branchMutations)):
                            if branchMutations[mutation2][2] == (branchMutations[mutation1][2] + 1):
                                positionsToRemove.append(mutation1)
                                positionsToRemove.append(mutation2)
                    if len(positionsToRemove) != 0:
                        #Write the positions to the mutations not analysed file
                        for removePosition in positionsToRemove:
                            outMutationsNotUsed.write(branchMutations[removePosition][0] + str(branchMutations[removePosition][1]) + branchMutations[removePosition][3] + "," + branchName + ",Double_substitution\n")
                        #Remove the positions from the mutations
                        for ele in sorted(positionsToRemove, reverse = True):
                            del branchMutations[ele]

                #Update the reference sequence to get the current context
                updatedReference = updateReference(tree, clade, branchMutationDict, referenceSequence)

                for mutation in branchMutations:
                    #Check if the mutation is at the end of the genome, if it is it does not have context so cannt be analysed
                    if (mutation[2] == 1) or (mutation[2] == referenceLength):
                        outMutationsNotUsed.write(mutation[0] + str(mutation[1]) + mutation[3] + "," + branchName + ",End_of_genome\n")
                    #Check if the mutation does not involve 2 nucleotides
                    elif (mutation[0] not in nucleotides) or (mutation[3] not in nucleotides):
                        outMutationsNotUsed.write(mutation[0] + str(mutation[1]) + mutation[3] + "," + branchName + ",Mutation_does_not_involve_two_nucleotides\n")
                    else:
                        mutationContext = getContext(mutation, updatedReference)

                        #Check if the upstream or downstream nucleotides are not A, C, G or T
                        if (mutationContext[0] not in nucleotides) or (mutationContext[1] not in nucleotides):
                            outMutationsNotUsed.write(mutation[0] + str(mutation[1]) + mutation[3] + "," + branchName + ",Surrounding_position_not_nucleotide\n")
                        else:
                            if (mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]) in spectrum:
                                spectrum[mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]] += 1
                            else:
                                spectrum[complement(mutationContext[0]) + complement(mutation[0]) + complement(mutation[3]) + complement(mutationContext[1])] += 1
    
    #Write the mutational spectrum
    for eachSpectrum in spectrum:
        outMutations.write(eachSpectrum + "," + str(spectrum[eachSpectrum]) + "\n")

    outMutations.close()
    outMutationsNotUsed.close()