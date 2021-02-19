#Combines multiple spectra from MutTui into a single catalog that can be submitted to signal for spectrum decomposition
#Can combine full mutational spectra or mutational type counts
#Each column in the combined catalog is a single input spectrum
#Can take an optional csv conversion file to rename samples in the output. This has 2 columns - column 1 is
#the path to the file, column 2 is its conversion name
#Also outputs a file containing the column names in the combined catalog and the corresponding original spectrum

import argparse
from plot_spectrum import convertSpectrumDict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        required = True,
                        nargs = "+",
                        help = "Set of spectra to be combined into a single catalog",
                        type = argparse.FileType("r"))
    parser.add_argument("-c",
                        "--conversion",
                        dest = "conversion",
                        required = False,
                        help = "Optional csv file containing sample name conversion. Should have " + 
                        "two columns with no header. Column 1 is the path to the file as it will " + 
                        "be given to the script. Column 2 is the name of the sample in the output",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outFile",
                        dest = "outFile",
                        required = True,
                        help = "Output csv catalog")
    
    args = parser.parse_args()

    outFile = open(args.outFile, "w")

    #List of spectra
    spectraDict = {}

    #Extract the spectra into spectraDict
    for spectrum in args.spectra:
        spectraDict[spectrum.name] = convertSpectrumDict(spectrum)
    
    #Dictionary with sample names as keys and files as values
    sampleDict = {}
    #If conversion provided, use this conversion
    if args.conversion:
        cFile = open(args.conversion.name).readlines()
        for line in cFile:
            sampleDict[line.strip().split(",")[1]] = line.strip().split(",")[0]
    #Otherwise name samples as Sample i
    else:
        for i, spectrum in enumerate(args.spectra):
            sampleDict["Sample" + str(i + 1)] = spectrum.name
    
    #Write the sample to file conversion
    conversionFile = open("sample_in_catalog_to_file_conversion.csv", "w")
    conversionFile.write("File,Sample_in_combined_catalog\n")
    for eachSample in sampleDict:
        conversionFile.write(sampleDict[eachSample] + "," + eachSample + "\n")
    
    #Write the header to the output file
    outFile.write("Substitution")
    for sample in sampleDict:
        outFile.write("," + sample)
    outFile.write("\n")
    
    #Iterate through the mutations and write their number in each sample
    firstKey = list(spectraDict.keys())[0]
    for mutation in spectraDict[firstKey]:
        outFile.write(mutation)
        for sample in sampleDict:
            outFile.write("," + str(spectraDict[sampleDict[sample]][mutation]))
        outFile.write("\n")

    outFile.close()
    conversionFile.close()