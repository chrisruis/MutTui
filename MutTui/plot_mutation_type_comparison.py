#Plots the proportion of each mutation type across a set of samples
#Runs on a multi-sample catalog which can be produced from a set of MutTui output spectra using combine_spectra_into_catalog.py
#This plots the 6 mutation types, not spectra

import argparse
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import sys

#Verifies that the mutations are in the expected order
def verifyOrder(catalog):
    inOrder = True

    if catalog[1].strip().split(",")[0] != "CA":
        inOrder = False
    if catalog[2].strip().split(",")[0] != "CG":
        inOrder = False
    if catalog[3].strip().split(",")[0] != "CT":
        inOrder = False
    if catalog[4].strip().split(",")[0] != "TA":
        inOrder = False
    if catalog[5].strip().split(",")[0] != "TC":
        inOrder = False
    if catalog[6].strip().split(",")[0] != "TG":
        inOrder = False
    
    return(inOrder)

#Converts the mutations in a catalog to proportions
def convertCatalogProportion(catalog):
    #Samples as keys, mutation type proportions as values
    sM = dict()

    for s in range(1, len(catalog[0].strip().split(","))):
        sampleName = catalog[0].strip().split(",")[s]
        tM = float(0)

        for line in catalog[1:]:
            tM += float(line.strip().split(",")[s])
        
        sM[sampleName] = list()

        for line in catalog[1:]:
            sM[sampleName].append(float(line.strip().split(",")[s])/tM)

    return(sM)

#Extracts colours from a colour file to a list
def getColours(colourFile, samples):
    colours = open(colourFile).readlines()

    colourDict = dict()
    for c in colours:
        colourDict[c.strip().split(",")[0]] = c.strip().split(",")[1]

    colourList = list()
    for s in colourDict:
        colourList.append(colourDict[s])
    
    return(colourList)

#Plots and saves a bar plot of mutation proportions
def plotProportions(xList, yList, colours, outFile):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.bar(xList, yList, color = colours)
    plt.xticks(rotation = "vertical")
    plt.subplots_adjust(bottom = 0.15)
    
    if type(outFile) == str:
        plt.savefig(outFile)
    else:
        plt.savefig(outFile.name)

#Plots the proportion of each mutation type in each of a set of samples
#The mutations are in order CA, CG, CT, TA, TC, TG, as verified by verifyOrder
def plotTypeComparison(catalogDict, colourFile):
    #Sample names for the plot
    samples = list()
    for s in catalogDict:
        samples.append(s)

    #Create lists for the proportion of each mutation type
    ca = list()
    cg = list()
    ct = list()
    ta = list()
    tc = list()
    tg = list()

    #Iterate through the samples and add to the lists
    for s in catalogDict:
        ca.append(catalogDict[s][0])
        cg.append(catalogDict[s][1])
        ct.append(catalogDict[s][2])
        ta.append(catalogDict[s][3])
        tc.append(catalogDict[s][4])
        tg.append(catalogDict[s][5])
    
    if not colourFile:
        plotProportions(samples, ca, "blue", "proportion_C_A.pdf")
        plotProportions(samples, cg, "black", "proportion_C_G.pdf")
        plotProportions(samples, ct, "red", "proportion_C_T.pdf")
        plotProportions(samples, ta, "grey", "proportion_T_A.pdf")
        plotProportions(samples, tc, "green", "proportion_T_C.pdf")
        plotProportions(samples, tg, "pink", "proportion_T_G.pdf")
    else:
        colours = getColours(colourFile, samples)
        plotProportions(samples, ca, colours, "proportion_C_A.pdf")
        plotProportions(samples, cg, colours, "proportion_C_G.pdf")
        plotProportions(samples, ct, colours, "proportion_C_T.pdf")
        plotProportions(samples, ta, colours, "proportion_T_A.pdf")
        plotProportions(samples, tc, colours, "proportion_T_C.pdf")
        plotProportions(samples, tg, colours, "proportion_T_G.pdf")

if __name__ == "__main__":
    description = "Plots mutation type proportions from a catalogue"
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument("-c",
                        "--catalog",
                        dest = "catalog",
                        required = True,
                        help = "Mutational catalog containing samples to be plotted. " +
                        "This should contain the six mutational types, not whole spectra, " +
                        "and can be generated using combine_spectra_into_catalog.py",
                        type = argparse.FileType("r"))
    parser.add_argument("-l",
                        "--colours",
                        dest = "colours",
                        default = False,
                        help = "Optional csv containing colours. Two columns with no header. " + 
                        "Column 1 is the sample name as it appears in the catalog, column 2 is the colour")
    
    args = parser.parse_args()

    #Import the catalog
    catalog = open(args.catalog.name).readlines()

    if verifyOrder(catalog) == False:
        sys.exit("Catalog is not in the correct order. Check the mutational order is CA, CG, CT, TA, TC, TG and try again")

    catalogProportions = convertCatalogProportion(catalog)

    plotTypeComparison(catalogProportions, args.colours)