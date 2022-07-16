#Functions to combine mutational spectra
#Typically run after MutTui to combine multiple calculated spectra into a single spectrum

import argparse
from .plot_spectrum import convertSpectrumDict

def combine_spectra_parser(parser):

    parser.description = "Sums the mutations in a number of input spectra to form a single output spectrum"

    parser.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        required = True,
                        nargs = "+",
                        help = "Spectrum files to be combined. All files specified with -s will be combined",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outFile",
                        dest = "outFile",
                        required = True,
                        help = "Name of the output file to which the combined spectrum will be written")

    parser.set_defaults(func=combine_spectra)

    return(parser)

def combine_spectra(args):

    #Import the first spectrum, convert it to a dictionary that the others will be combined into
    combinedDict = convertSpectrumDict(args.spectra[0])

    outFile = open(args.outFile, "w")
    outFile.write("Substitution,Number_of_mutations\n")

    #Iterate through the remaining spectra and add them to combinedDict
    for spectrumFile in args.spectra[1:]:
        spectrum = convertSpectrumDict(spectrumFile)

        for mutation in combinedDict:
            combinedDict[mutation] += spectrum[mutation]
    
    #Write the combined spectrum
    for eachMutation in combinedDict:
        outFile.write(eachMutation + "," + str(combinedDict[eachMutation]) + "\n")

    outFile.close()

    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = combine_spectra_parser(parser)
    args = parser.parse_args()

    # run combine_spectra
    args.func(args)

    return

if __name__ == "__main__":
    main()
