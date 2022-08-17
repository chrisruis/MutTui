#Functions for comparing spectra
#Run on multiple spectra output from MutTui and/or combined with combine_spectra.py
#Converts each spectrum from mutation counts to mutation proportions and compares the proportions
#Can carry out statistical tests to identify significantly different mutations between 2 spectra

import argparse
import numpy as np
from collections import Counter
from reconstruct_spectrum import getMutationDict
from plot_spectrum import *

#from .reconstruct_spectrum import getMutationDict
#from .plot_spectrum import *

#Imports the mutation counts from a spectrum file
def importMC(sFile):
    #All mutations
    s = list()

    #Unique mutations
    mSet = set()

    with open(sFile.name) as f:
        next(f)
        for l in f:
            s += [l.strip().split(",")[0].split("[")[1].split("]")[0]] * int(l.strip().split(",")[1])
            mSet.add(l.strip().split(",")[0])
    
    s = np.array(s, dtype = "<U3")
    
    return(s, mSet)

#Permutes the values within 2 lists
def permuteValues(allP, numMut1):
    #Shuffle the mutations
    rs = np.random.choice(allP, len(allP), replace = False)
    #rs = random.sample(allP.keys(), tm, counts = allP.values())

    #Assign mutations to samples
    rs1 = Counter(rs[:int(numMut1)])
    rs2 = Counter(rs[int(numMut1):len(rs)])

    return(rs1, rs2)

#Compares mutation type proportions between spectra using a permutation test
def compareMutationTypes(spectrum1, spectrum2, rna, out):
    #Import spectra
    s1, mSet1 = importMC(spectrum1)
    s2, mSet2 = importMC(spectrum2)

    #Verify the mutations are the same in each spectrum
    if mSet1 != mSet2:
        raise RuntimeError("Spectra do not contain the same mutations")

    #Combine mutations
    allP = np.concatenate((s1, s2))
    
    outF = open(out + "_mutation_type_comparison.csv", "w")
    outF.write("Mutation,Difference,p-value\n")

    #Mutation types, only first 6 used for DNA spectra
    mutTypes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "G>T", "G>C", "G>A", "A>T", "A>G", "A>C"]
    #Number of mutations to iterate through
    if rna:
        nMT = 12
    else:
        nMT = 6
    
    #Number of mutations in each spectrum
    numMut1 = len(s1)
    numMut2 = len(s2)
    tm = numMut1 + numMut2

    #Differences in real data
    sDiff = dict()
    #Lists of differences in permutations
    pDiff = dict()

    #Count occurrences of each mutation type in input spectra
    c1 = Counter(s1)
    c2 = Counter(s2)
    #Iterate through the mutations, add to sDiff and pDiff and calculate difference in real data
    for m in mutTypes[:nMT]:
        m1 = float(c1[m])/numMut1
        m2 = float(c2[m])/numMut2
        sDiff[m] = m1 - m2

        pDiff[m] = list()

    #Iterate through permutations, for each permute values, iterate through the mutations,
    #calculate their difference in the permutation and count permutations with at least the same difference
    for eP in range(1000):
        #Permute mutations
        rs1, rs2 = permuteValues(allP, numMut1)

        #Calculate difference for each mutation type and add to pDiff
        for m in mutTypes[:nMT]:
            pM1 = float(rs1[m])/numMut1
            pM2 = float(rs2[m])/numMut2
            pDiff[m].append(pM1 - pM2)

    #Count permutations with at least the same difference
    for m in mutTypes[:nMT]:
        numP = 0
        for i in pDiff[m]:
            if i >= sDiff[m]:
                numP += 1

        outF.write(m + "," + str(sDiff[m]) + "," + str((numP + 1)/1001) + "\n")
    
    outF.close()

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
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen",
                        action = "store_true",
                        default = False)
    parser.add_argument("--compare_mutations",
                        dest = "cm",
                        help = "Compare the mutation type proportions between spectra. This runs a permutation test " + 
                        "on each mutation type to compare the proportion difference between spectra with that expected by chance. An " + 
                        "additional output file out_prefix_mutation_type_comparison.csv is written. Should be used without --proportions",
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
    
    #Compare mutation type proportions if specified
    if args.cm:
        compareMutationTypes(args.spectra[0], args.spectra[1], args.rna, args.out_prefix)

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

    if not args.rna:
        plotSpectrumComparison(spectrumDifference, args.out_prefix + "_spectrum_comparison.pdf")
        plotSpectrumPointComparison(spectrum1Proportions, spectrum2Proportions, args.out_prefix + "_proportion_comparison.pdf")
    else:
        plotRNASpectrumComparison(spectrumDifference, args.out_prefix + "_spectrum_comparison.pdf")

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