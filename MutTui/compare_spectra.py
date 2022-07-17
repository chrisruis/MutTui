#Functions for comparing spectra
#Run on multiple spectra output from MutTui and/or combined with combine_spectra.py
#Converts each spectrum from mutation counts to mutation proportions and compares the proportions
#Can carry out statistical tests to identify significantly different mutations between 2 spectra

import argparse
from .reconstruct_spectrum import getMutationDict
from .plot_spectrum import *

#Converts a given spectrum from number of mutations to proportion of mutations and returns the spectrum as a dictionary
#def convertSpectrumProportions(spectrum):
#    totalMutations = float(0)
#
#    for mutation in spectrum:
#        totalMutations += float(spectrum[mutation])
#    
#    for eachMutation in spectrum:
#        spectrum[eachMutation] = float(spectrum[eachMutation])/totalMutations
#       
#    return(spectrum)

#Converts a given spectrum from number of mutations to proportion of mutations and returns a list of proportions
#def convertSpectrumProportions(spectrum):
#    totalMutations = float(0)
#
#    for mutation in spectrum:
#        totalMutations += float(spectrum[mutation])
#    
#    sL = list()
#    
#    for eachMutation in spectrum:
#        sL.append(float(spectrum[eachMutation])/totalMutations)
#    
#    return(sL)

def compare_spectra_parser(parser):

    parser.description = "Subtracts one mutational spectrum from another"

    parser.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        required = True,
                        nargs = 2,
                        help = "The 2 spectra that will be compared. The second provided spectrum will be " + 
                        "subtracted from the first. Provide both spectra to -s separated by a space",
                        type = argparse.FileType("r"))
    parser.add_argument("--proportions",
                        dest = "proportions",
                        help = "Specify if the spectra to be compared are proportions rather than numbers of mutations, " + 
                        "e.g. if comparing SigProfilerExtractor signatures",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "--out",
                        dest = "out_prefix",
                        required = True,
                        help = "The prefix of the output files")

    parser.set_defaults(func=compare_spectra)

    return(parser)

def compare_spectra(args):

    #Import the spectra into dictionaries
    if args.proportions:
        spectrum1Proportions = convertSpectrumDictProportions(args.spectra[0])
        spectrum2Proportions = convertSpectrumDictProportions(args.spectra[1])
    else:
        spectrum1 = convertSpectrumDict(args.spectra[0])
        spectrum2 = convertSpectrumDict(args.spectra[1])
        #Convert the spectra to proportions
        spectrum1Proportions = convertSpectrumProportions(spectrum1)
        spectrum2Proportions = convertSpectrumProportions(spectrum2)

    #Will be filled with the difference between each mutation type
    spectrumDifference = {}

    #Iterate through the mutations and add to spectrumDifference
    for mutation in spectrum1Proportions:
        spectrumDifference[mutation] = spectrum1Proportions[mutation] - spectrum2Proportions[mutation]
    
    outFile = open(args.out_prefix + ".csv", "w")
    outFile.write("Substitution,Difference_in_spectra\n")

    for eachMutation in spectrumDifference:
        outFile.write(eachMutation + "," + str(spectrumDifference[eachMutation]) + "\n")

    outFile.close()

    plotSpectrumComparison(spectrumDifference, args.out_prefix + "_spectrum_comparison.pdf")

    plotSpectrumPointComparison(spectrum1Proportions, spectrum2Proportions, args.out_prefix + "_proportion_comparison.pdf")

    return

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = compare_spectra_parser(parser)
    args = parser.parse_args()

    # run compare_spectra
    args.func(args)

    return

if  __name__ == "__main__":
    main()