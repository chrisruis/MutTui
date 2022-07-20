#Calculates the distance between all pairs of mutational spectra and clusters based on these distances
#Calculates distance using either the Bhattacharyya distance or Jensen-Shannon distance
#Can optionally take a file with information to colour the points
#The colouring file has 2 columns - column 1 is the file paths as they will be given to the script, column 2 is the colours

import os
import argparse
import numpy as np
from scipy import spatial
from sklearn import manifold, decomposition
import pandas as pd
import matplotlib
matplotlib.use('AGG')
from matplotlib import pyplot as plt
from .plot_spectrum import convertSpectrumDictProportions

#Converts a set of spectra into a list of dictionaries
def convertSpectraList(spectra, labels):
    #List of spectra dictionaries
    spectraList = [{} for i in range(len(spectra))]

    sampleNames = list()

    #Iterate through the spectra, convert to proportions and add to spectraList
    for i, spectrum in enumerate(spectra):
        if labels:
            sampleNames.append(labels[i])
        else:
            sampleNames.append(spectrum.name)
        
        #Calculate the total number of mutations in the spectrum
        totalMutations = float(0)
        mD = convertSpectrumDictProportions(spectrum)
        for m in mD:
            totalMutations += mD[m]
        
        #Calculate the proportion of each mutation and add to spectraList and sL
        for m in mD:
            spectraList[i][m] = mD[m]/totalMutations
    
    return(spectraList, sampleNames)

#Converts a multi-sample catalog into a list of dictionaries
def convertCatalog(catalogFile):
    catalog = open(catalogFile.name).readlines()

    sampleNames = catalog[0].strip().split(",")[1:]

    #List of spectra dictionaries
    spectraList = [{} for i in range(len(sampleNames))]

    #Calculate the total mutations in each sample if the spectra are not proportions
    totalMutations = [0] * len(sampleNames)
    for row in catalog[1:]:
        for sample in range(len(sampleNames)):
            totalMutations[sample] += float(row.strip().split(",")[sample + 1])
    
    #Calculate the proportion of each mutation in each sample
    for row in catalog[1:]:
        for sample in range(len(sampleNames)):
            spectraList[sample][row.strip().split(",")[0]] = float(row.strip().split(",")[sample + 1])/totalMutations[sample]
    
    return(spectraList, sampleNames)

#Converts a given colouring file to a dictionary with files as keys and colours as values and
#an array of colours for umap
def getColourDict(colourFile):
    if colourFile is None:
        return(None, None)
    
    else:
        colourDict = {}

        #Conversion from label to identifier
        conversion = dict()
        #Incremented with each new label
        iterator = 0
        colourList = []

        #Iterate through the labels, add the label to colourList
        #Add to the label conversion if not already present
        with open(colourFile.name) as fileobject:
            for line in fileobject:
                l = line.strip().split("\t")[1]
                colourDict[line.strip().split("\t")[0]] = l
                if l not in conversion:
                    conversion[l] = iterator
                    iterator += 1
                colourList.append(conversion[l])
        
        return(colourDict, colourList)

#Calculates the cosine distance between 2 given spectra
#Used to allow multiple distance metrics, now only supports cosine distance
def calculateSpectraDistance(spectrum1, spectrum2):
    s1 = list(spectrum1.values())
    s2 = list(spectrum2.values())
    distance = spatial.distance.cosine(s1, s2)
    
    return(distance)

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
            distances[spectrum1, spectrum2] = calculateSpectraDistance(spectraList[spectrum1], spectraList[spectrum2])
            distances[spectrum2, spectrum1] = distances[spectrum1, spectrum2]
    
    return(distances)

#Converts a list of spectrum dictionaries to an array
def arraySpectra(spectra):
    #Extract mutations from keys
    m = spectra[0].keys()

    sa = np.zeros((len(spectra), 96))
    
    for i, s in enumerate(spectra):
        for j, em in enumerate(m):
            sa[i,j] = s[em]
    
    return(sa)

#PCA of SBS spectra based on mutation proportions
def sbsPCA(spectra, sn, col, cl, output_dir):
    #PCA of spectra
    pca = decomposition.PCA(n_components = 2)
    pc = pca.fit_transform(spectra)

    #Extract coordinates to plot
    x = list()
    y = list()
    for i in pc:
        x.append(i[0])
        y.append(i[1])
    
    #Colours for plot
    if col is None:
        colours = ["black" for i in range(len(sn))]
    else:
        colours = []
        for s in sn:
            colours.append(col[s])
    
    #Write PCA coordinates
    with open(output_dir + "SBS_PCA_coordinates.csv", "w") as outP:
        outP.write("Sample,Principle_component_1,Principle_component_2\n")
        for es in range(len(sn)):
            outP.write(sn[es] + "," + str(x[es]) + "," + str(y[es]) + "\n")
    
    #Data frame for plotting
    pDF = pd.DataFrame({"x": x, "y": y, "cols": colours})
    
    #Plot PCA
    plt.style.use("ggplot")
    fig = plt.figure()
    #Colour column is labels
    if col and not cl:
        groups = pDF.groupby("cols")
        for name, group in groups:
            plt.scatter(group.x, group.y, label = name)
    #Colour column is colours
    else:
        plt.scatter(pDF.x, pDF.y, color = pDF.cols)
    plt.grid(True)
    plt.xlabel("Principle component 1")
    plt.ylabel("Principle component 2")
    plt.tight_layout()
    fig.savefig(output_dir + "SBS_PCA.pdf")

#Plots a multidimensional scaling based on a given distance matrix
def plotMDS(distances, file_names, colourDict, output_dir):
    #MDS of the distances
    mds = manifold.MDS(n_components = 2, dissimilarity = "precomputed")
    projection = mds.fit(distances)
    coords = projection.embedding_

    #Extract the colours of each point
    #colourDict = getColourDict(colourFile)
    colours = []
    if colourDict is None:
        for eachFile in file_names:
            colours.append("blue")
    else:
        for eachFile in file_names:
            #Removed .name from after eachFile as it seems like this will have strings as keys
            #if a colour file is provided
            colours.append(colourDict[eachFile])

    #Write MDS coordinates
    with open(output_dir + "distance_mds_coordinates.txt", "w") as points_out:
        points_out.write("Sample\tx_coordinate\ty_coordinate\n")
        for i, coord in zip(file_names, coords):
            points_out.write("%s\t%s\t%s\n" % (i, coord[0], coord[1]))
    
    #Identify plot margins
    c_min = np.min(coords) - abs(np.quantile(coords, 0.05))
    c_max = np.max(coords) + abs(np.quantile(coords, 0.05))

    #Plot MDS
    plt.style.use("ggplot")
    fig = plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1], color = colours)
    plt.grid(True)
    plt.xlabel("MDS dimension 1")
    plt.ylabel("MDS dimension 2")
    plt.xlim((c_min, c_max))
    plt.ylim((c_min, c_max))
    plt.tight_layout()
    fig.savefig(output_dir + "distance_MDS.pdf")

#OLD FUNCTION USED TO PLOT UMAP BASED ON SAMPLE DISTANCES
#UMAP NO LONGER PLOTTED
#UMAP clustering and plotting
def plotUMAP(distances, colours, output_dir):
    distanceMap = umap.UMAP(n_neighbors = 5, metric = "cosine").fit(distances)

    if colours is not None:
        labels = np.array(colours)
        umapFig = umap.plot.points(distanceMap, labels = labels)
    else:
        umapFig = umap.plot.points(distanceMap)

    umapFig.figure.savefig(output_dir + "sample_umap.pdf")

#OLD FUNCTION USED TO RUN UMAP ON A SET OF SPECTRA
#UMAP NO LONGER USED
#Calculates a UMAP embedding for a set of spectra
def extractUMAP(sL, labels, colours, output_dir):
    #Convert to array
    sL = np.asarray(sL)

    reducer = umap.UMAP()

    embedding = reducer.fit_transform(sL)

    #Write the sample coordinates
    out_coordinates = open(output_dir + "sample_umap_coordinates.csv", "w")
    out_coordinates.write("Sample,x_coordinate,y_coordinate\n")
    for l in range(len(embedding)):
        out_coordinates.write(labels[l] + "," + str(embedding[l][0]) + "," + str(embedding[l][1]) + "\n")
    out_coordinates.close()

    fig = plt.figure()
    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.gca().set_aspect("equal", "datalim")
    fig.savefig(output_dir + "test.pdf")

    distanceMap = umap.UMAP().fit(sL)
    if colours is not None:
        labels = np.array(colours)
        umapFig = umap.plot.points(distanceMap, labels = labels)
    else:
        umapFig = umap.plot.points(distanceMap)
    
    umapFig.figure.savefig(output_dir + "sample_umap.pdf")


def cluster_spectra(args):
    #If using labels, verify that there are the same number of labels as spectra
    if args.labels:
        if len(args.spectra) != len(args.labels):
            raise RuntimeError("The number of labels provided with -l must match the number of spectra provided with -s")
    
    #Make sure trailing forward slash is present in output directory
    args.output_dir = os.path.join(args.output_dir, "")

    #If multiple spectra are provided, extract these to lists
    if args.spectra:
        spectraList, sampleNames = convertSpectraList(args.spectra, args.labels)
    else:
        #Extract the catalog into spectraList
        spectraList, sampleNames = convertCatalog(args.catalog)
    
    #Calculate distances between all pairs of spectra
    distances = getDistanceMatrix(spectraList)

    #Write cosine similarity and distances
    s_out = open(args.output_dir + "cosine_similarity.csv", "w")
    d_out = open(args.output_dir + "cosine_distances.csv", "w")
    s_out.write("Sample")
    d_out.write("Sample")
    for sample in sampleNames:
        s_out.write("," + sample)
        d_out.write("," + sample)
    s_out.write("\n")
    d_out.write("\n")
    for row in range(len(distances)):
        s_out.write(sampleNames[row])
        d_out.write(sampleNames[row])
        for column in distances[row]:
            s_out.write("," + str(1 - column))
            d_out.write("," + str(column))
        s_out.write("\n")
        d_out.write("\n")
    s_out.close()
    d_out.close()
    
    #Extract the colours from the colour file if present
    colourDict, cConversion = getColourDict(args.colour_file)

    #Extract spectra to arrays
    sA = arraySpectra(spectraList)

    #MDS of mutation proportions
    sbsPCA(sA, sampleNames, colourDict, args.colour_labels, args.output_dir)

    #MDS of sample distances
    if args.mds_distance:
        plotMDS(distances, sampleNames, colourDict, args.output_dir)

    return

def cluster_spectra_parser(parser):

    parser.description = "Calculate distances and cluster a set of input SBS spectra. These can be provided as individual spectrum files with -s or as a combined catalogue with -c"

    spectra = parser.add_mutually_exclusive_group(required = True)
    spectra.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        nargs = "+",
                        help = "SBS spectrum files to be clustered. All files specified with -s will be clustered",
                        type = argparse.FileType("r"))
    spectra.add_argument("-c",
                        "--catalog",
                        dest = "catalog",
                        help = "Multi-sample catalogue containing SBS spectra, samples as columns and mutation counts/proportions as rows",
                        type = argparse.FileType("r"))
    
    parser.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        nargs = "+",
                        help = "Labels for each spectrum provided with -s. If not provided, the file names will be used as labels",
                        default = None)
    parser.add_argument("-cl",
                        "--colours",
                        dest = "colour_file",
                        default = None,
                        help = "Optional file containing information to colour points in the output clustering. " + 
                        "This file should contain 2 columns separated by tabs with no header. Column 1 is the " + 
                        "name of the sample. If providing labels with -l, these names should match the labels. If " +
                        "using a catalog, these names should match the column names. If using file paths, these names should " +
                        "match the file paths, as provided to -s. Column 2 contains information for colouring. This can either be " + 
                        "the group to which the sample name belongs (e.g. enteric or group1) or a colour for each sample. If " + 
                        "using colours, add the --colour_labels option",
                        type = argparse.FileType("r"))
    parser.add_argument("--mutation_types",
                        dest = "mt",
                        help = "Run PCA on mutation proportions within each of the 6 mutation types. This will output " + 
                        "6 additional csv files and 6 additional plots, one per mutation type in each case. This calculates " + 
                        "the proportion of each contextual mutation within the respective mutation type (e.g. the proportion of " + 
                        "C>A mutations that are T[C>A]G) ands runs a PCA on those 16 proportions across spectra",
                        action = "store_true",
                        default = False)
    parser.add_argument("--plot_mds_distance",
                        dest = "mds_distance",
                        help = "Output an MDS plot based on sample distances",
                        action = "store_true",
                        default = False)
    parser.add_argument("--colour_labels",
                        dest = "colour_labels",
                        help = "Specify that the sample labels in the file provided with -cl are colours",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "--out_prefix",
                        dest = "output_dir",
                        required = True,
                        help = "Output directory")

    parser.set_defaults(func=cluster_spectra)

    return(parser)



def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = cluster_spectra_parser(parser)
    args = parser.parse_args()

    # run cluster_spectra
    args.func(args)

    return

if  __name__ == "__main__":
    main()