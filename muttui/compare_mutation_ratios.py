#Compares the ratios of each mutation type between a set of mutational spectra
#Run on mutation_types_label_X.csv files from MutTui
#Outputs the ratio between the number of each pair of mutation types in each spectrum

import argparse
from plot_spectrum import convertSpectrumDict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--spectra",
                        dest = "spectra",
                        required = True,
                        nargs = "+",
                        help = "mutation_types_label.csv files from MutTui to be compared",
                        type = argparse.FileType("r"))
    parser.add_argument("-l",
                        "--labels",
                        dest = "labels",
                        nargs = "+",
                        help = "Label of each spectrum, one per spectrum. Samples will be labelled with this " +
                        "in the output file. Default is to name each spectrum as it is given to -s",
                        default = None)
    parser.add_argument("-o",
                        help = "Name of output file")
    args = parser.parse_args()

    outFile = open(args.o, "w")
    outFile.write("Spectrum,C>A/C>G,C>A/C>T,C>A/T>A,C>A/T>C,C>A/T>G,C>G/C>T,C>G/T>A,C>G/T>C,C>G/T>G,C>T/T>A,C>T/T>C,C>T/T>G,T>A/T>C,T>A/T>G,T>C/T>G\n")

    #Set the sample labels
    if args.labels:
        if len(args.spectra) != len(args.labels):
            raise RuntimeError("The number of labels provided with -l needs to match the number of spectra provided with -s. Please check the labels")
        labels = args.labels
    else:
        labels = list()
        for s in args.spectra:
            labels.append(s.name)
    
    #Iterate through the spectra, calculate the ratios of each mutation type and write
    for i, s in enumerate(args.spectra):
        spectrum = convertSpectrumDict(s)
        outFile.write(labels[i] + "," + str(spectrum["CA"]/spectrum["CG"]) + "," + str(spectrum["CA"]/spectrum["CT"]) + "," + str(spectrum["CA"]/spectrum["TA"]) + "," + str(spectrum["CA"]/spectrum["TC"]) + "," + str(spectrum["CA"]/spectrum["TG"]) + "," + str(spectrum["CG"]/spectrum["CT"]) + "," + str(spectrum["CG"]/spectrum["TA"]) + "," + str(spectrum["CG"]/spectrum["TC"]) + "," + str(spectrum["CG"]/spectrum["TG"]) + "," + str(spectrum["CT"]/spectrum["TA"]) + "," + str(spectrum["CT"]/spectrum["TC"]) + "," + str(spectrum["CT"]/spectrum["TG"]) + "," + str(spectrum["TA"]/spectrum["TC"]) + "," + str(spectrum["TA"]/spectrum["TG"]) + "," + str(spectrum["TC"]/spectrum["TG"]) + "\n")

    outFile.close()