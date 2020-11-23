#Calculates the distance between all pairs of mutational spectra and clusters based on these distances
#Calculates distance using either the Bhattacharyya distance or Jensen-Shannon distance
#Can optionally take a file with information to colour the points
#The colouring file has 2 columns - column 1 is the file paths as they will be given to the script, column 2 is the colours

import os
import argparse
from math import sqrt, log2
import numpy as np
from sklearn import manifold
from matplotlib import pyplot as plt
from plot_spectrum import convertSpectrumDict
from compare_spectra import convertSpectrumProportions

#Converts a given colouring file to a dictionary with files as keys and colours as distances
def getColourDict(colourFile, spectra):
    print(colourFile)

#Calculates the distance between 2 given spectra
def calculateSpectraDistance(spectrum1, spectrum2, method):
    distance = float(0)

    if method == "Bhattacharyya":
        for mutation in spectrum1:
            distance += sqrt(spectrum1[mutation] * spectrum2[mutation])
        Bhattacharyya_distance = -log2(distance)
    
    return(Bhattacharyya_distance)

#Creates a matrix of zeros with number of rows and columns equal to the provided number of samples
def getZerosMatrix(nSamples):
    dist_mat = np.zeros((nSamples, nSamples), dtype = float)

    return(dist_mat)

#Calculates a distance matrix between all pairs of spectra
def getDistanceMatrix(spectraList):
    #Distance matrix of zeros
    distances = getZerosMatrix(len(spectraList))

    #Iterate through the pairs of spectra and calculate their distance
    for spectrum1 in range(len(spectraList) - 1):
        for spectrum2 in range((spectrum1 + 1), len(spectraList)):
            distances[spectrum1, spectrum2] = calculateSpectraDistance(spectraList[spectrum1], spectraList[spectrum2], args.method)
            distances[spectrum2, spectrum1] = distances[spectrum1, spectrum2]
    
    return(distances)

#Plots a multidimensional scaling based on a given distance matrix
def plotMDS(distances, file_names, output_dir):
    #MDS of the distances
    mds = manifold.MDS(n_components = 2, dissimilarity = "precomputed")
    projection = mds.fit(distances)
    coords = projection.embedding_

    #Write MDS coordinates
    with open(output_dir + "mds_coordinates.txt", "w") as points_out:
        points_out.write("Sample\tx_coordinate\ty_coordinate\n")
        for i, coord in zip(file_names, coords):
            points_out.write("%s\t%s\t%s\n" % (i.name, coord[0], coord[1]))
    
    #Identify plot margins
    c_min = np.min(coords) - abs(np.quantile(coords, 0.05))
    c_max = np.max(coords) + abs(np.quantile(coords, 0.05))

    #Plot MDS
    plt.style.use("ggplot")
    fig = plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1])
    plt.grid(True)
    plt.xlabel("MDS dimension 1")
    plt.ylabel("MDS dimension 2")
    plt.xlim((c_min, c_max))
    plt.ylim((c_min, c_max))
    plt.tight_layout()
    fig.savefig(output_dir + "sample_MDS.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        required = True,
                        nargs = "+",
                        help = "Spectrum files to be clustered. All files specified with -s will be clustered",
                        type = argparse.FileType("r"))
    parser.add_argument("-m",
                        "--method",
                        dest = "method",
                        default = "Bhattacharyya",
                        help = "The method used to calculate distances between pairs of spectra. " + 
                        "Options are Bhattacharyya (default) and JS (Jensen-Shannon)")
    parser.add_argument("-c",
                        "--colours",
                        dest = "colour_file",
                        default = None,
                        help = "Optional file containing information to colour points in the output clustering. " + 
                        "This file should contain 2 columns separated by tabs with no header. Column 1 is the " + 
                        "file path of each spectrum that will be clustered. Give the path from the directory in " + 
                        "which you will run the script as it would be provided to -s. Column 2 is the colour that the " + 
                        "corresponding point will be.")
    parser.add_argument("-o",
                        "--out_prefix",
                        dest = "output_dir",
                        required = True,
                        help = "Output directory")
    
    args = parser.parse_args()

    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #Extract the colours to a dictionary if required, otherwise set the points to blue
    getColourDict(args.colourFile, args.spectra)

    #List of spectra
    spectraList = []

    #Extract the spectra into spectraDict
    for spectrum in args.spectra:
        spectraList.append(convertSpectrumProportions(convertSpectrumDict(spectrum)))
    
    #Calculate distances between all pairs of spectra
    distances = getDistanceMatrix(spectraList)

    #Write distances
    with open(args.output_dir + "sample_distances.txt", "w") as distances_out:
        distances_out.write("Sample")
        for sample in args.spectra:
            distances_out.write("\t" + sample.name)
        distances_out.write("\n")
        for row in range(len(distances)):
            distances_out.write(args.spectra[row].name)
            for column in distances[row]:
                distances_out.write("\t" + str(column))
            distances_out.write("\n")

    plotMDS(distances, args.spectra, args.output_dir)