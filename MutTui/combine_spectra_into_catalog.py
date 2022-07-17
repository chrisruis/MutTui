#Combines multiple spectra from MutTui into a single catalog that can be submitted to signal for spectrum decomposition
#Can combine full mutational spectra or mutational type counts
#Each column in the combined catalog is a single input spectrum
#Can take an optional csv conversion file to rename samples in the output. This has 2 columns - column 1 is
#the path to the file, column 2 is its conversion name
#Also outputs a file containing the column names in the combined catalog and the corresponding original spectrum

import argparse
from .plot_spectrum import convertSpectrumDict, convertSpectrumDictProportions


def combine_spectra_into_catalog(args):

    #Check a maximum of 1 label option is provided
    if args.conversion and args.labels:
        raise RuntimeError("Please provide either labels with -l or a label file with -c, not both")

    outFile = open(args.outFile + ".csv", "w")
    outFileMatrix = open(args.outFile + ".txt", "w")

    #List of spectra
    spectraDict = {}

    #Extract the spectra into spectraDict
    for spectrum in args.spectra:
        if args.proportions:
            spectraDict[spectrum.name] = convertSpectrumDictProportions(spectrum)
        else:
            spectraDict[spectrum.name] = convertSpectrumDict(spectrum)
    
    #Dictionary with sample names as keys and files as values
    sampleDict = {}
    #If conversion provided, use this conversion
    if args.conversion:
        cFile = open(args.conversion.name).readlines()
        for line in cFile:
            sampleDict[line.strip().split("\t")[1]] = line.strip().split("\t")[0]
    #If labels are provided, use these
    elif args.labels:
        if len(args.spectra) != len(args.labels):
            raise RuntimeError("The number of labels needs to match the number of spectra")
        for i in range(len(args.spectra)):
            sampleDict[args.labels[i]] = args.spectra[i].name
    #Otherwise name samples as Sample i
    else:
        for i, spectrum in enumerate(args.spectra):
            sampleDict["Sample" + str(i + 1)] = spectrum.name
    
    #Write the sample to file conversion
    conversionFile = open("sample_in_catalog_to_file_conversion.csv", "w")
    conversionFile.write("File,Sample_in_combined_catalog\n")
    for eachSample in sampleDict:
        conversionFile.write(sampleDict[eachSample] + "," + eachSample + "\n")
    
    #Write the header to the output files
    outFile.write("Substitution")
    outFileMatrix.write("Mutation Types")
    for sample in sampleDict:
        outFile.write("," + sample)
        outFileMatrix.write("\t" + sample)
    outFile.write("\n")
    outFileMatrix.write("\n")
    
    #Iterate through the mutations and write their number in each sample
    firstKey = list(spectraDict.keys())[0]
    for mutation in spectraDict[firstKey]:
        outFile.write(mutation)
        outFileMatrix.write(mutation)
        for sample in sampleDict:
            outFile.write("," + str(spectraDict[sampleDict[sample]][mutation]))
            outFileMatrix.write("\t" + str(spectraDict[sampleDict[sample]][mutation]))
        outFile.write("\n")
        outFileMatrix.write("\n")

    outFile.close()
    conversionFile.close()

    return


def combine_spectra_into_catalog_parser(parser):

    parser.description = "Combines multiple input spectra and/or catalogues into a single output catalogue"


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
                        help = "Optional tab separated file containing sample name conversion. Should have " + 
                        "two columns with no header. Column 1 is the path to the file as it will " + 
                        "be given to the script. Column 2 is the name of the sample in the output",
                        type = argparse.FileType("r"),
                        default = False)
    parser.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        nargs = "+",
                        help = "Optional labels to be given to spectra in the combined catalogue. Needs to have " + 
                        "one label for each spectrum provided with -s",
                        default = False)
    parser.add_argument("--proportions",
                        dest = "proportions",
                        help = "Specify if the spectra are proportions rather than numbers of mutations, " + 
                        "e.g. if using SigProfilerExtractor signatures",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "--outFile",
                        dest = "outFile",
                        required = True,
                        help = "Output file prefix. Two files will be written - a csv catalog" + 
                        " and a tab separated mutation matrix")
    
    parser.set_defaults(func=combine_spectra_into_catalog)

    return(parser)



def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = combine_spectra_into_catalog_parser(parser)
    args = parser.parse_args()

    # run combine_spectra_into_catalog
    args.func(args)

    return

if __name__ == "__main__":
    main()
