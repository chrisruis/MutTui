#Functions to plot mutational spectra
#Used by muttui.py to plot spectra but can also be run as a standalone script on a MutTui output spectrum

import argparse
import matplotlib.pyplot as plt

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
    
    plt.style.use("ggplot")
    plt.bar(mutations, mutationCounts, color = colours)
    plt.xlabel("Mutation")
    plt.ylabel("Number of mutations")
    plt.tick_params(axis = "x", which = "both", bottom = False, labelbottom = False)
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
    
    args = parser.parse_args()

    #Extract the spectrum to a dictionary with mutations as keys and counts as values
    spectrum = convertSpectrumDict(args.spectrum_file)

    plotSpectrumFromDict(spectrum, args.outFile)