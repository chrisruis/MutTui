#Splits the mutations in an all_included_mutations.csv file from MutTui into branch groups
#Branch groups are defined with a branch labelling csv file with 2 columns, first column is the
#branch names as they appear in the tree, second column is the label of the branch

import argparse
import pandas as pd
import os
from isvalid import *
from reconstruct_spectrum import *
from plot_spectrum import *

#Extracts the labels to a dictionary
def getLabels(labelsFile):
    labels = dict()

    with open(labelsFile.name) as f:
        for l in f:
            labels[l.strip().split(",")[0]] = l.strip().split(",")[1]
    
    return(labels)

#Post-processes mutations into branch groups
def postProcessMutations(mutationsFile, labelsFile, rna, outdir):
    #Import branch labels and extract to dictionary
    labels = getLabels(labelsFile)

    #Generate an empty spectrum and mutation dictionary for each label
    spectraDict = dict()
    allMDict = dict()
    for label in set(labels.values()):
        if not rna:
            spectraDict[label] = getMutationDict()
        else:
            spectraDict[label] = getRNADict()
        allMDict[label] = list()
    
    #Iterate through the mutations and add to the respective spectrum
    mutations = pd.read_csv(mutationsFile.name)
    for i in range(mutations.shape[0]):
        branch = mutations["Branch"][i]
        if branch in labels:
            bL = labels[branch]
            spectraDict[bL][mutations["Substitution"][i].replace("[", "").replace(">", "").replace("]", "")] += 1
            allMDict[bL].append([mutations["Mutation_in_alignment"][i], mutations["Mutation_in_genome"][i], mutations["Substitution"][i], mutations["Branch"][i]])
    
    #Write and plot the spectra
    for eachLabel in spectraDict:
        outFile = open(outdir + "mutational_spectrum_label_" + eachLabel + ".csv", "w")
        outFile.write("Substitution,Number_of_mutations\n")
        for eachMutation in spectraDict[eachLabel]:
            outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(spectraDict[eachLabel][eachMutation]) + "\n")
        outFile.close()

        #Plot the spectrum
        outSpectrum = open(outdir + "mutational_spectrum_label_" + eachLabel + ".pdf", "w")
        spectrumFormat = convertSpectrumFormat(spectraDict[eachLabel])
        if not rna:
            plotSpectrumFromDict(spectrumFormat, False, outSpectrum)
        else:
            plotRNA(spectrumFormat, False, outSpectrum)
        outSpectrum.close()
    
        #Write the mutations in each label
        outMutation = open(outdir + "all_included_mutations_label_" + eachLabel + ".csv", "w")
        outMutation.write("Mutation_in_alignment,Mutation_in_genome,Substitution,Branch\n")
        for eM in allMDict[eachLabel]:
            outMutation.write(",".join(eM) + "\n")
        outMutation.close()

#Filters mutations with fewer than a given number of occurrences and within a given region of the genome
def filterMutations(mutationsFile, nm, genome_pos, align_pos, rna, outdir):
    #Set max and min genome and alignment positions
    if not genome_pos:
        g_pos = [1, 100000000]
    else:
        g_pos = [int(genome_pos[0]), int(genome_pos[1])]
    if not align_pos:
        a_pos = [1, 100000000]
    else:
        a_pos = [int(align_pos[0]), int(align_pos[1])]
    
    #Max mutations to keep
    mMut = int(nm)

    #Count occurrences of each mutation
    mutations = dict()
    with open(mutationsFile.name) as f:
        next(f)
        for l in f:
            #Extract genome and alignment positions
            g = int(l.strip().split(",")[1][1:-1])
            a = int(l.strip().split(",")[0][1:-1])

            #Count mutations within required genome and alignment regions
            if (g >= g_pos[0]) and (g <= g_pos[1]) and (a >= a_pos[0]) and (a <= a_pos[1]):
                if l.strip().split(",")[0] not in mutations:
                    mutations[l.strip().split(",")[0]] = 0
                mutations[l.strip().split(",")[0]] += 1
    
    #Empty spectrum
    if rna:
        spectrum = getRNADict()
    else:
        spectrum = getMutationDict()
    
    #Filter mutations
    with open(mutationsFile.name) as f:
        next(f)
        for l in f:
            if l.strip().split(",")[0] in mutations:
                if mutations[l.strip().split(",")[0]] <= mMut:
                    substitution = l.strip().split(",")[2]
                    spectrum[substitution[0] + substitution[2] + substitution[4] + substitution[6]] += 1
    
    #Write spectrum
    outFile = open(outdir + "mutational_spectrum_filtered.csv", "w")
    outFile.write("Substitution,Number_of_mutations\n")
    for eachMutation in spectrum:
        outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(spectrum[eachMutation]) + "\n")
    outFile.close()

    #Plot spectrum
    outSpectrum = open(outdir + "mutational_spectrum_filtered.pdf", "w")
    spectrumFormat = convertSpectrumFormat(spectrum)
    if not rna:
        plotSpectrumFromDict(spectrumFormat, False, outSpectrum)
    else:
        plotRNA(spectrumFormat, False, outSpectrum)
    outSpectrum.close()

if __name__ == "__main__":
    description = "Splits mutations into separate mutational spectra based on branch labels"
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument("-m",
                        "--mutations",
                        dest = "mutations",
                        help = "all_included_mutations.csv file from MutTui",
                        type = argparse.FileType("r"))
    parser.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        help = "csv file used to split mutations into categories. This should have 2 columns: column 1 " + 
                        "is the branch name as it appears in all_included_mutations.csv and column 2 is its label. This file should not " + 
                        "have a header",
                        type = argparse.FileType("r"))
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen",
                        action = "store_true",
                        default = False)
    parser.add_argument("--filter",
                        dest = "filter",
                        help = "Specify to filter mutations based on their number of occurrences. If specified, use -n " + 
                        "to specify the maximum number of occurrences to keep a mutation",
                        action = "store_true",
                        default = False)
    parser.add_argument("-n",
                        "--number_mutations",
                        dest = "number_mutations",
                        help = "The maximum number of occurrences to keep a mutation. Mutation inferred to occur on more occasions than this " + 
                        "will be filtered out. Default 1000000 to not exclude mutations by occurrences. Only used when --filter is specified",
                        default = "1000000")
    parser.add_argument("-g",
                        "--genome_positions",
                        dest = "genome_positions",
                        help = "The first and last positions in the genome to be included when filtering. Mutations outside of these genome " + 
                        "positions will be excluded. Defaults set 1 and 100,000,000 to include all positions. Only used when --filter is specified",
                        nargs = 2,
                        default = None)
    parser.add_argument("-a",
                        "--align_positions",
                        dest = "align_positions",
                        help = "The first and last positions in the alignment to be included when filtering. Mutations outside of these alignment " + 
                        "positions will be excluded. Defaults set 1 and 100,000,000 to include all positions. Only used when --filter is specified",
                        nargs = 2,
                        default = None)
    parser.add_argument("-o",
                        "--out_dir",
                        dest = "output_dir",
                        help = "Location of output directory, should already be created and ideally be empty",
                        type = lambda x: is_valid_folder(parser, x))
    args = parser.parse_args()

    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #Post-process the spectrum into groups
    if not args.filter:
        postProcessMutations(args.mutations, args.labels, args.rna, args.output_dir)
    #Filter the mutations
    else:
        filterMutations(args.mutations, args.number_mutations, args.genome_positions, args.align_positions, args.rna, args.output_dir)