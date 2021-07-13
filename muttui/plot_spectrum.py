#Functions to plot mutational spectra
#Used by muttui.py to plot spectra but can also be run as a standalone script on a MutTui output spectrum

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#Converts an input spectrum from a csv file to a dictionary
def convertSpectrumDict(spectrumFile):
    spectrum = open(spectrumFile.name).readlines()

    spectrumDict = {}

    for mutation in spectrum[1:]:
        spectrumDict[mutation.strip().split(",")[0]] = int(mutation.strip().split(",")[1])
    
    return(spectrumDict)

#Converts a spectrum to the required format from plotting
#Converts from format ACAA to A[C>A]A
def convertSpectrumFormat(spectrum):
    spectrumDict = {}

    for mutation in spectrum:
        spectrumDict[mutation[0] + "[" + mutation[1] + ">" + mutation[2] + "]" + mutation[3]] = spectrum[mutation]
    
    return(spectrumDict)

#Plots a spectrum from a dictionary with mutations as keys and mutation counts as values
def plotSpectrumFromDict(spectrum, outFile):
    #The 96 DNA mutations as a list, used so they are always in the same order
    mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    #Will be filled with the mutation counts
    mutationCounts = []

    #Iterate through the mutations and add their count to mutationCounts
    for mutation in mutations:
        mutationCounts.append(spectrum[mutation])
    
    #Colours of the bars
    colourSet = ["blue", "black", "red", "grey", "green", "pink"]
    colours = [i for i in colourSet for j in range(16)]

    #Labels of the mutation type rectangles
    labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

    #Used to plot the rectangles above the plot containing the mutation type
    rect_lower = float(np.max(mutationCounts)) + float(np.max(mutationCounts) * 0.05)
    rect_width = float(np.max(mutationCounts)) * 0.1

    #The coordinates of the mutation type rectangles and text
    mutation_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    rect_coords = [-0.5, 15.5, 31.5, 47.5, 63.5, 79.5]
    text_coords = [4, 20, 36, 52, 68, 84]

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.bar(mutations, mutationCounts, color = colours)
    #ax.add_patch(plt.Rectangle((80, rect_lower), 15, rect_width, clip_on = False))
    for i, rect in enumerate(rect_coords):
        ax.add_patch(plt.Rectangle((rect, rect_lower), 16, rect_width, facecolor = colourSet[i]))
        ax.text(text_coords[i], (rect_lower + (rect_width/3)), mutation_types[i], color = "white", fontweight = "bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels([""] * len(mutations))
    plt.xlabel("Mutation")
    plt.ylabel("Number of mutations")
    plt.margins(0)
    if type(outFile) == str:
        plt.savefig(outFile)
    else:
        plt.savefig(outFile.name)

#Plots the number of mutations of each type, combining all contexts from that type
def plotMutationType(mtCounts, outFile):
    #Extract the mutations and counts to separate lists
    mutations = list(mtCounts.keys())
    mutationCounts = list(mtCounts.values())
    #Convert the counts to proportions
    #mutationProportions = list()
    #Used to plot proportion of mutations, now changed to plotting number of mutations so not used
    #totalMutations = float(sum(mutationCounts))
    #for m in mutationCounts:
    #    mutationProportions.append(float(m)/totalMutations)

    #Colours of the bars
    colourSet = ["blue", "black", "red", "grey", "green", "pink"]

    fig = plt.figure()
    ax = plt.subplot(111)
    #Used to plot proportion of mutations, now changed to plotting number of mutations so not used
    #ax.bar(mutations, mutationProportions, color = colourSet)
    ax.bar(mutations, mutationCounts, color = colourSet)
    plt.xlabel("Mutation type")
    plt.ylabel("Proportion of mutations")

    if type(outFile) == str:
        plt.savefig(outFile)
    else:
        plt.savefig(outFile.name)

#Plots a comparison of the proportions of each mutation in 2 spectra
#x-axis is the proportion of the mutation in sample 1, y-axis is the proportion of the mutation in sample 2
#Takes 2 spectra to be compared and the output PDF to save to
#The input spectra should have already been converted to proportions
def plotSpectrumPointComparison(spectrum1, spectrum2, outFile):
    #The 96 DNA mutations as a list, used so they are always in the same order
    mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    #Colours of the points
    colourSet = ["blue", "black", "red", "grey", "green", "pink"]
    colours = [i for i in colourSet for j in range(16)]

    #The mutation proportions in spectrum 1
    mutation1 = []
    #The mutation proportions in spectrum 2
    mutation2 = []

    for mutation in mutations:
        mutation1.append(spectrum1[mutation])
        mutation2.append(spectrum2[mutation])
    
    fig = plt.figure()
    ax = plt.subplot(111)
    scatter = ax.scatter(mutation1, mutation2, color = colours, s = 8)
    plt.xlabel("Mutation proportion in sample 1")
    plt.ylabel("Mutation proportion in sample 2")
    plt.margins(0)

    #Calculate plot limits and plot x=y line across the plot
    lims = [0, np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, color = "black", alpha = 0.75, linestyle = "--")
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    #Add legend
    mutation_legend = [Line2D([0], [0], color = "blue", lw = 2),
                    Line2D([0], [0], color = "black", lw = 2),
                    Line2D([0], [0], color = "red", lw = 2),
                    Line2D([0], [0], color = "grey", lw = 2),
                    Line2D([0], [0], color = "green", lw = 2),
                    Line2D([0], [0], color = "pink", lw = 2)]
    ax.legend(mutation_legend, ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"], bbox_to_anchor = (0.8, 1.15), ncol = 3)

    #ax.legend(, loc = "upper right")

    #plt.legend(bbox_to_anchor = (1.05, 1))

    if type(outFile == str):
        plt.savefig(outFile)
    else:
        plt.savefig(outFile.name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s",
                        "--spectrum",
                        dest = "spectrum_file",
                        required = True,
                        help = "csv file containing spectrum from MutTui",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outFile",
                        dest = "outFile",
                        required = True,
                        help = "Output PDF file to which spectrum will be written")
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen, will plot an RNA mutational spectrum",
                        action = "store_true",
                        default = False)
    
    args = parser.parse_args()

    #Extract the spectrum to a dictionary with mutations as keys and counts as values
    spectrum = convertSpectrumDict(args.spectrum_file)

    plotSpectrumFromDict(spectrum, args.outFile)