#Calculates the distance between all pairs of mutational spectra and clusters based on these distances
#Calculates distance using either the Bhattacharyya distance or Jensen-Shannon distance

import argparse
from math import sqrt, log2
from plot_spectrum import convertSpectrumDict
from compare_spectra import convertSpectrumProportions

#Calculates the distance between 2 given spectra
def calculateSpectraDistance(spectrum1, spectrum2, method):
    distance = float(0)

    if method == "Bhattacharyya":
        for mutation in spectrum1:
            distance += sqrt(spectrum1[mutation] * spectrum2[mutation])
        Bhattacharyya_distance = -log2(distance)
    
    return(Bhattacharyya_distance)

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
    parser.add_argument("-o",
                        "--out_prefix",
                        required = True,
                        help = "Output file prefix")
    
    args = parser.parse_args()

    #List of spectra
    spectraList = []

    #Extract the spectra into spectraDict
    for spectrum in args.spectra:
        spectraList.append(convertSpectrumProportions(convertSpectrumDict(spectrum)))
    
    #Iterate through the pairs of spectra and calculate their distance
    for spectrum1 in range(len(spectraList) - 1):
        for spectrum2 in range((spectrum1 + 1), len(spectraList)):
            print(calculateSpectraDistance(spectraList[spectrum1], spectraList[spectrum2], args.method))