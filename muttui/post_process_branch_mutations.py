#Splits the mutations in an all_included_mutations.csv file from MutTui into branch groups
#Branch groups are defined with a branch labelling csv file with 2 columns, first column is the
#branch names as they appear in the tree, second column is the label of the branch

import argparse
from label_tree_gene_mutations import labelTreeMutations
import pandas as pd
import os
from isvalid import *
from reconstruct_spectrum import *
from plot_spectrum import *

#Extracts the labels to a dictionary
def getLabels(labelsFile):
    labels = dict()

    l = pd.read_csv(labelsFile.name)

    for i in range(l.shape[0]):
        labels[l["branch"][i]] = l["label"][i]
    
    return(labels)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m",
                        "--mutations",
                        dest = "mutations",
                        required = True,
                        help = "all_included_mutations.csv file from MutTui",
                        type = argparse.FileType("r"))
    parser.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        required = True,
                        help = "csv file used to split mutations into categories. This should have 2 columns: column 1 " + 
                        "is the branch name as it appears in all_included_mutations.csv and column 2 is its label",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--out_dir",
                        dest = "output_dir",
                        required = True,
                        help = "Location of output directory, should already be created and ideally be empty",
                        type = lambda x: is_valid_folder(parser, x))
    args = parser.parse_args()

    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #Import branch labels and extract to dictionary
    labels = getLabels(args.labels)

    #Generate an empty spectrum for each label
    spectraDict = dict()
    for label in set(labels.values()):
        spectraDict[label] = getMutationDict()
    
    #Iterate through the mutations and add to the respective spectrum
    mutations = pd.read_csv(args.mutations.name)
    for i in range(mutations.shape[0]):
        spectraDict[labels[mutations["Branch"][i]]][mutations["Substitution"][i].replace("[", "").replace(">", "").replace("]", "")] += 1
    
    #Write and plot the spectra
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