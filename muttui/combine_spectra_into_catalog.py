#Combines multiple spectra from MutTui into a single catalog that can be submitted to signal for spectrum decomposition
#Each column in the combined catalog is a single input spectrum
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
    for mutation in spectraDict[sampleDict["Sample1"]]:
        outFile.write(mutation)
        for sample in sampleDict:
            outFile.write("," + str(spectraDict[sampleDict[sample]][mutation]))
        outFile.write("\n")

    outFile.close()
    conversionFile.close()