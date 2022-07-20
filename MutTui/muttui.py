#Runs the MutTui pipeline on an input tree and alignment

import os
import argparse
from .treetime import run_treetime, change_gaps_to_Ns
from .isvalid import *
import subprocess
import array
from Bio import AlignIO, Phylo
from collections import OrderedDict
import gffutils as gff
from .add_tree_node_labels import *
from .branch_labelling import *
from .reconstruct_spectrum import *
from .plot_spectrum import *
from .gff_conversion import *

from .__init__ import __version__




def muttui_parser(parser):

    parser.description = "Run the MutTui pipeline on a given alignment and tree"

    #Options for input and output files
    io_opts = parser.add_argument_group("Input/output")
    io_opts.add_argument("-a",
                        "--alignment",
                        dest = "alignment",
                        required = True,
                        help = "Input fasta alignment",
                        type = argparse.FileType("r"))
    io_opts.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick tree file",
                        type = argparse.FileType("r"))
    io_opts.add_argument("-o",
                        "--out_dir",
                        dest = "output_dir",
                        required = True,
                        help = "Location of output directory, should already be created and ideally be empty",
                        type = lambda x: is_valid_folder(parser, x))
    io_opts.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        help = "Comma separated file with sequences and their clade labels to divide the spectrum. " +
                        "If this option is included, treetime mugration will be run to identify changes in the label across the tree. " +
                        "Does not need to include all taxa. Taxa not included will be given a label of OTHER by default, " +
                        "or OTHER_A if OTHER is already taken. This file should have a header.",
                        default = None)
    io_opts.add_argument("-lt",
                        "--labelled_tree",
                        dest = "labelled_tree",
                        help = "Tree with state labels, should be the same tree as that provided with -t but with nodes " + 
                        "labelled with their state from label_tree.py",
                        default = None,
                        type = argparse.FileType("r"))
    io_opts.add_argument("-r",
                        "--reference",
                        dest = "reference",
                        help = "Reference genome sequence containing all sites, used to get context of mutations, not required if using --all_sites",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-c",
                        "--conversion",
                        dest = "conversion",
                        help = "Conversion file for alignment position to genome position, used to get context of mutations, not required if using --all_sites",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-g",
                        "--gff",
                        dest = "gff",
                        help = "GFF reference containing gene coordinates in reference sequence. Used to split mutations into transcription strands " + 
                        "and identify synonymous mutations when --synonymous is used",
                        type = argparse.FileType("r"),
                        default = None)
    io_opts.add_argument("-to",
                        "--treetime_out",
                        dest = "treetime_out",
                        help = "The location of the directory containing treetime output files from ancestral reconstruction. Only " + 
                        "used if option --start_from_treetime is specified, in which case the output files in this directory will " + 
                        "be used to calculate the mutational spectrum",
                        default = None)
    io_opts.add_argument("-bm",
                        "--branch_mutations",
                        dest = "branch_mutations",
                        help = "The minimum number of mutations on a branch to calculate its mutational spectrum. Only used when " + 
                        "specifying --branch_specific which will output the spectrum for each branch with more than this number of " + 
                        "mutations separately. Default = 50",
                        default = "50")
    
    #Options to include in treetime reconstruction
    treetime = parser.add_argument_group("treetime command")
    treetime.add_argument("--add_treetime_cmds",
                        dest = "add_treetime_cmds",
                        help = "Additional options to supply to treetime (these are not checked). Supply these together in quotes",
                        type = str,
                        default = None)
    
    #Other options
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen, MutTui will output an RNA mutational spectrum",
                        action = "store_true",
                        default = False)
    parser.add_argument("--all_sites",
                        dest = "all_sites",
                        help = "Specify that the alignment contains all sites, in which case a reference genome does not need to be provided",
                        action = "store_true",
                        default = False)
    parser.add_argument("--include_all_branches",
                        dest = "include_all_branches",
                        help = "Use when specifying a labelled tree with -lt. By default, the branches along which the label " + 
                        "changes are excluded as it is often not clear at what point along the branch the label changed. With " + 
                        "adding --include_all_branches, all of the branches along which a label changes will be included in the " + 
                        "spectrum of the downstream label",
                        action = "store_true",
                        default = False)
    parser.add_argument("--exclude_root_branches",
                        dest = "exclude_root_branches",
                        help = "Use when excluding the two branches that branch directly from the root node. These branches will " +
                        "not be included in any spectrum",
                        action = "store_true",
                        default = False)
    parser.add_argument("-rs",
                        "--root_state",
                        dest = "root_state",
                        help = "Specify the root state of the given label if known. " + 
                        "The root state is calculated by treetime mugration and in almost all cases will be resolved. " + 
                        "If the root state has multiple assignments that are equally likely, it cannot be assigned from " + 
                        "treetime mugration output. In these rare cases, use this option to assign the root state. This option " + 
                        "can also be used to assign a root state if you'd like but its recommended to use the mugration state",
                        default = None)
    parser.add_argument("--start_from_treetime",
                        dest = "start_from_treetime",
                        help = "Use this option to start with treetime output and so skip inference of ancestral mutations. Use this " + 
                        "if you have already run treetime. The directory containing the treetime output files needs to be provided with -to",
                        action = "store_true",
                        default = False)
    parser.add_argument("--strand_bias",
                        dest = "strand_bias",
                        help = "Split mutations into transcribed and untranscribed strands to test for transcription strand bias. A GFF " + 
                        "file will need to be provided with option -g to identify genes",
                        action = "store_true",
                        default = False)
    parser.add_argument("--synonymous",
                        dest = "synonymous",
                        help = "Use non-coding and synonymous mutations only to calculate the mutational spectrum. A GFF file will need to be " + 
                        "provided with option -g which will be used to identify genes",
                        action = "store_true",
                        default = False)
    parser.add_argument("--branch_specific",
                        dest = "branch_specific",
                        help = "Calculate the mutational spectrum for each branch in the tree separately. Only branches containing at least " + 
                        "the number of mutations specified with -bm will be included",
                        action = "store_true",
                        default = False)
    parser.add_argument("--version",
                        action = "version",
                        version = "%(prog)s " + __version__)

    parser.set_defaults(func=muttui)

    return parser



def muttui(args):

    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #If --branch_specific is used, set --include_all_branches to true
    if args.branch_specific:
        args.include_all_branches = True
    
    #Minimum mutations to include a branch, only used with --branch_specific
    minimumMutations = int(args.branch_mutations)

    #Open output files
    outMutationsNotUsed = open(args.output_dir + "mutations_not_included.csv", "w")
    outMutationsNotUsed.write("Mutation_in_alignment,Mutation_in_genome,Branch,Reason_not_included\n")

    outAllMutations = open(args.output_dir + "all_included_mutations.csv", "w")
    outAllMutations.write("Mutation_in_alignment,Mutation_in_genome,Substitution,Branch\n")

    outAllDouble = open(args.output_dir + "all_included_double_substitutions.csv", "w")
    outAllDouble.write("Mutation_in_alignment,Mutation_in_genome,Substitution,Branch,Original_mutation\n")

    if not args.start_from_treetime:
        print("Running treetime ancestral reconstruction to identify mutations")

        #Convert gaps to Ns in the alignment, run treetime on the new alignment
        change_gaps_to_Ns(args.alignment, args.output_dir)
        run_treetime(open(args.output_dir + "gaps_to_N_alignment.fasta"), args.tree, args.output_dir, args.add_treetime_cmds)
    
        print("treetime reconstruction complete. Importing alignment from reconstruction and tree")

        #Import the alignment from treetime
        alignment = AlignIO.read(args.output_dir + "ancestral_sequences.fasta", "fasta")
        
    else:
        #Import the already calculated alignment from treetime
        print("Importing alignment from reconstruction and tree")
        #Make sure the trailing forward slash is present in output directory
        args.treetime_out = os.path.join(args.treetime_out, "")

        #Import the alignment from treetime
        alignment = AlignIO.read(args.treetime_out + "ancestral_sequences.fasta", "fasta")
    
    #Import the original unlabelled tree
    tree = Phylo.read(args.tree.name, "newick")
    #Ladderize the tree so the branches are in the same order as the treetime tree
    tree.ladderize()
    tree = labelBranchesTreetime(tree)

    print("Alignment and tree imported. Reconstructing spectrum")

    #Check if a GFF file is needed, if so read it in and process it
    if args.strand_bias or args.synonymous:
        if not args.gff:
            raise RuntimeError("GFF file needs to be provided with -g when using --strand_bias or --synonymous")
        else:
            geneCoordinates, positionGene = convertGFF(args.gff.name)

    #Label branches in the tree into categories, each category will have a separate spectrum
    if args.labels:
        ####The code in this if statement has not been altered with branch_mutations.txt
        ####The elif and else sections have been altered
        #Extract the labels to a dictionary and add taxa without a label
        labelDict = getLabelDict(tree, args.labels)

        #Write the labels to a csv that can be used by treetime mugration
        writeLabels(labelDict, args.output_dir)

        #Run treetime mugration to reconstruct the clade labels across the tree
        run_treetime_mugration(args.output_dir + "annotated_tree.nexus", args.output_dir + "all_taxon_labels.csv", args.output_dir)

        #Import the files from treetime mugration
        mugrationTree = Phylo.read(args.output_dir + "mugration_out/annotated_tree.nexus", "nexus")
        confidence = open(args.output_dir + "mugration_out/confidence.csv").readlines()
        gtr = open(args.output_dir + "mugration_out/GTR.txt").readlines()

        #Identify the root state and add it to the mugrationTree
        mugrationTree = rootState(mugrationTree, confidence, gtr, args.root_state)

        labelledTree, treeLabels = labelBranchesMugration(tree, mugrationTree)
    #Label branches from the provided tree
    elif args.labelled_tree:
        labelledTree, treeLabels = getLabelledTreeLabels(tree, args.labelled_tree)
    #Label branches with their branch name
    elif args.branch_specific:
        labelledTree, treeLabels = labelBranchesNames(tree)
    #Label all branches with the same label
    else:
        labelledTree, treeLabels = labelAllBranches(tree)
    
    #Branch categories as keys, spectra as values
    spectraDict = {}
    doubleSpectraDict = {}
    if args.strand_bias:
        sbDict = {}
    #Create empty spectrum for each branch category
    if args.rna:
        for label in treeLabels:
            spectraDict[label] = getRNADict()
    else:
        for label in treeLabels:
            spectraDict[label] = getMutationDict()
            doubleSpectraDict[label] = getDoubleSubstitutionDict()
            if args.strand_bias:
                sbDict[label] = getStrandBiasDict()
    
    #The 4 nucleotides, used to check if mutated, upstream and downstream bases are nucleotides
    nucleotides = ["A","C","G","T"]

    #Convert the positions in the alignment to genome positions, if --all_sites specified the positions will be the same
    if args.all_sites:
        positionTranslation = allSitesTranslation(alignment)
    else:
        positionTranslation = convertTranslation(args.conversion)
    
    #Extracts mutations to a dictionary from the branch_mutations.txt file
    if not args.start_from_treetime:
        branchMutationDict = getBranchMutationNexusDict(args.output_dir + "annotated_tree.nexus", positionTranslation)
    else:
        branchMutationDict = getBranchMutationNexusDict(args.output_dir + "annotated_tree.nexus", positionTranslation)

    #Get the reference sequence, if -r specified this will be the provided genome, otherwise all sites in the alignment are assumed
    #and the root sequence from the ancestral reconstruction is used
    referenceSequence = getReference(args.reference, args.all_sites, alignment, positionTranslation)
    referenceLength = len(referenceSequence)
    
    #Iterate through the branches, get the category of the branch, identify the contextual mutations, add to the corresponding spectrum
    for clade in labelledTree.find_clades():
        #Check if there are mutations along the current branch, only need to analyse branches with mutations
        if clade.name in branchMutationDict:
            #The label of the current branch, this will be None if the label changes along this branch
            branchCategory = getBranchCategory(labelledTree, clade, args.include_all_branches)

            #Extract the mutations along the branch. This will be None if there are no mutations but treetime has still added the mutation comment
            branchMutations = branchMutationDict[clade.name]

            #If using --branch_specific, check if the branch contains at least -bm mutations, if so
            #add the branch to spectraDict. If not, set branchCategory to None so it won't be analysed
            if args.branch_specific:
                if len(branchMutations) >= minimumMutations:
                    if args.rna:
                        spectraDict[clade.name] = getRNADict()
                    else:
                        spectraDict[clade.name] = getMutationDict()
                else:
                    branchCategory = None
            
            #If using --exclude_root_branches, check if the branch comes off the root, if so set the branchCategory
            #to None so it won't be analysed
            if args.exclude_root_branches:
                if len(labelledTree.get_path(clade)) == 1:
                    branchCategory = None

            #Check if the branch has a category, will be None if the branch is a transition between categories
            #and option --include_all_branches is not specified
            if branchCategory is not None:
                #Extract double substitutions, remove mutations at the ends of the genome or not involving 2 nucleotides
                branchMutations, doubleSubstitutions = filterMutations(branchMutations, clade, nucleotides, referenceLength, outMutationsNotUsed)
                
                #Update the reference sequence to get the current context
                updatedReference = updateReference(tree, clade, branchMutationDict, referenceSequence)

                #Check if only synonymous mutations should be included, if so filter the mutations
                if args.synonymous:
                    branchMutations = extractSynonymous(branchMutations, updatedReference, geneCoordinates, positionGene)
                
                for mutation in branchMutations:
                    mutationContext = getContext(mutation, updatedReference)
                    
                    #Check if the upstream or downstream nucleotides are not A, C, G or T
                    if (mutationContext[0] not in nucleotides) or (mutationContext[1] not in nucleotides):
                        outMutationsNotUsed.write(mutation[0] + str(mutation[1]) + mutation[3] + "," + mutation[0] + str(mutation[2]) + mutation[3] + "," + clade.name + ",Surrounding_position_not_nucleotide\n")
                    else:
                        #This will be true for all RNA mutations and half of DNA mutations
                        if (mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]) in spectraDict[branchCategory]:
                            spectraDict[branchCategory][mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]] += 1
                            outAllMutations.write(mutation[0] + str(mutation[1]) + mutation[3] + "," + mutation[0] + str(mutation[2]) + mutation[3] + "," + mutationContext[0] + "[" + mutation[0] + ">" + mutation[3] + "]" + mutationContext[1] + "," + clade.name + "\n")
                        #Add to the corresponding complement
                        else:
                            spectraDict[branchCategory][complement(mutationContext[1]) + complement(mutation[0]) + complement(mutation[3]) + complement(mutationContext[0])] += 1
                            outAllMutations.write(complement(mutation[0]) + str(mutation[1]) + complement(mutation[3]) + "," + complement(mutation[0]) + str(mutation[2]) + complement(mutation[3]) + "," + complement(mutationContext[1]) + "[" + complement(mutation[0]) + ">" + complement(mutation[3]) + "]" + complement(mutationContext[0]) + "," + clade.name + "\n")
                    
                        #If calculating strand bias, check if the mutation position is in a gene, if so extract its strand bias
                        if (args.strand_bias) and (mutation[2] in positionGene):
                            sb = getStrandBias(mutation, updatedReference, geneCoordinates, positionGene)
                            #Exclude genes in multiple genes on different strands
                            if sb is not None:
                                #Add to the corresponding dictionary
                                if (sb + mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]) in sbDict[branchCategory]:
                                    sbDict[branchCategory][sb + mutationContext[0] + mutation[0] + mutation[3] + mutationContext[1]] += 1
                                else:
                                    sbDict[branchCategory][sb + complement(mutationContext[1]) + complement(mutation[0]) + complement(mutation[3]) + complement(mutationContext[0])] += 1

                #Add double substitutions to the corresponding spectrum
                if len(doubleSubstitutions) > 0:
                    dsIter = iter(doubleSubstitutions)
                    for s1, s2 in zip(dsIter, dsIter):
                        #Check if both mutations involve nucleotides, if not write both to the unused file
                        if (s1[0] not in nucleotides) or (s1[3] not in nucleotides) or (s2[0] not in nucleotides) or (s2[3] not in nucleotides):
                            outMutationsNotUsed.write(s1[0] + str(s1[1]) + s1[3] + "," + s1[0] + str(s1[2]) + s1[3] + "," + clade.name + ",Double_substitution_does_not_involve_two_nucleotides\n")
                            outMutationsNotUsed.write(s2[0] + str(s2[1]) + s2[3] + "," + s2[0] + str(s2[2]) + s2[3] + "," + clade.name + ",Double_substitution_does_not_involve_two_nucleotides\n")
                        else:
                            #This will be true for all RNA mutations and half of DNA mutations
                            if (s1[0] + s2[0] + s1[3] + s2[3]) in doubleSpectraDict[branchCategory]:
                                doubleSpectraDict[branchCategory][s1[0] + s2[0] + s1[3] + s2[3]] += 1
                                outAllDouble.write(s1[0] + s2[0] + str(s1[1]) + s1[3] + s2[3] + "," + s1[0] + s2[0] + str(s1[2]) + s1[3] + s2[3] + "," + s1[0] + s2[0] + ">" + s1[3] + s2[3] + "," + clade.name + ",Forward\n")
                            else:
                                doubleSpectraDict[branchCategory][complement(s2[0]) + complement(s1[0]) + complement(s2[3]) + complement(s1[3])] += 1
                                outAllDouble.write(complement(s2[0]) + complement(s1[0]) + str(s1[1]) + complement(s2[3]) + complement(s1[3]) + "," + complement(s2[0]) + complement(s1[0]) + str(s1[2]) + complement(s2[3]) + complement(s1[3]) + "," + complement(s2[0]) + complement(s1[0]) + ">" + complement(s2[3]) + complement(s1[3]) + "," + clade.name + ",Reverse\n")
            #Write the mutations on the branch to the not used file
            else:
                for eM in branchMutations:
                    outMutationsNotUsed.write(eM[0] + str(eM[1]) + eM[3] + "," + eM[0] + str(eM[2]) + eM[3] + "," + clade.name + ",Label_changes_on_branch\n")
    
    #Write the spectra to separate files
    for eachLabel in spectraDict:
        outFile = open(args.output_dir + "mutational_spectrum_label_" + eachLabel + ".csv", "w")
        outFile.write("Substitution,Number_of_mutations\n")
        for eachMutation in spectraDict[eachLabel]:
            outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(spectraDict[eachLabel][eachMutation]) + "\n")
        outFile.close()

        #Plot the spectrum
        outSpectrum = open(args.output_dir + "mutational_spectrum_label_" + eachLabel + ".pdf", "w")
        spectrumFormat = convertSpectrumFormat(spectraDict[eachLabel])
        plotSpectrumFromDict(spectrumFormat, outSpectrum)
        outSpectrum.close()

        #Calculate the number of each type of mutation
        mtCounts = mutationTypeCount(spectraDict[eachLabel], args.rna)
        #Write the mutation type counts
        outMT = open(args.output_dir + "mutation_types_label_" + eachLabel + ".csv", "w")
        outMT.write("Mutation_type,Number_of_mutations\n")
        for m in mtCounts:
            outMT.write(m + "," + str(mtCounts[m]) + "\n")
        outMT.close()

        #Plot the mutation type counts
        outMTSpectrum = open(args.output_dir + "mutation_types_label_" + eachLabel + ".pdf", "w")
        plotMutationType(mtCounts, outMTSpectrum)
        outMTSpectrum.close()

        #Write the double substitution spectrum
        outDouble = open(args.output_dir + "DBS_label_" + eachLabel + ".csv", "w")
        outDouble.write("Substitution,Number_of_mutations\n")
        for eachMutation in doubleSpectraDict[eachLabel]:
            outDouble.write(eachMutation[:2] + ">" + eachMutation[2:] + "," + str(doubleSpectraDict[eachLabel][eachMutation]) + "\n")
        
        #Plot the double substitution spectrum
        outDoubleSpectrum = open(args.output_dir + "DBS_label_" + eachLabel + ".pdf", "w")
        doubleFormat = convertSpectrumFormat(doubleSpectraDict[eachLabel])
        plotDouble(doubleFormat, False, outDoubleSpectrum)
        outDoubleSpectrum.close()

        #Write the strand bias spectrum
        if args.strand_bias:
            outSB = open(args.output_dir + "strand_bias_label_" + eachLabel + ".csv", "w")
            outSB.write("Strand,Substitution,Number_of_mutations\n")
            for sbm in sbDict[eachLabel]:
                outSB.write(sbm[0] + "," + sbm[1:] + "," + str(sbDict[eachLabel][sbm]) + "\n")
            outSB.close()
        
            #Plot the strand bias spectrum
            outSBP = open(args.output_dir + "strand_bias_label_" + eachLabel + ".pdf", "w")
            plotSB(sbDict[eachLabel], False, outSBP)
            outSBP.close()
    
    #Write the spectra to a combined catalog if there is more than 1 label
    if len(spectraDict.keys()) > 1:
        outCombined = open(args.output_dir + "combined_catalog.csv", "w")
        outCombined.write("Substitution")
        for eachLabel in spectraDict:
            outCombined.write("," + eachLabel)
        outCombined.write("\n")

        for eachMutation in spectraDict[list(spectraDict.keys())[0]]:
            outCombined.write(eachMutation)
            for eachLabel in spectraDict:
                outCombined.write("," + str(spectraDict[eachLabel][eachMutation]))
            outCombined.write("\n")

        outCombined.close()

    #Close output files
    outMutationsNotUsed.close()
    outAllMutations.close()
    outAllDouble.close()

    return



def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = muttui_parser(parser)
    args = parser.parse_args()

    # run muttui
    args.func(args)

    return

if __name__ == "__main__":
    main()