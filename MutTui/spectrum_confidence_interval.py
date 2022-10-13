#Calculates confidence intervals for all mutations in a spectrum using Wilson score interval
#Can run on SBS, DBS or mutation type spectra

import argparse
from math import sqrt
from plot_spectrum import convertSpectrumDictProportions

#Calculates Wilson score interval for a mutation proportion
def wilsonSI(p, n):
    z = 1.96

    if (p == 0) or (n == 0):
        return(0, 0)

    denominator = 1 + z**2/n
    centre_adjusted_probability = p + z*z / (2*n)
    adjusted_standard_deviation = sqrt((p*(1 - p) + z*z / (4*n)) / n)
    lower_bound = (centre_adjusted_probability - z*adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z*adjusted_standard_deviation) / denominator

    return(lower_bound, upper_bound)

if __name__ == "__main__":
    description = "Calculates confidence intervals for all mutations in a SBS, DBS or mutation type spectrum"
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument("-s",
                        "--spectrum",
                        dest = "spectrum",
                        required = True,
                        help = "Mutational spectrum, can be SBS, DBS or mutation type",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Name of output csv file containing mutation proportions and confidence intervals")
    
    args = parser.parse_args()

    #Import the spectrum and convert to a dictionary
    spectrum = convertSpectrumDictProportions(args.spectrum)

    #Sum mutations
    tm = sum(spectrum.values())

    #Convert spectum to proportions
    sP = dict()
    for m in spectrum:
        sP[m] = spectrum[m]/tm
    
    out = open(args.outFile, "w")
    out.write("Substitution,Proportion_of_mutations,CI_lower,CI_upper\n")

    #Calculate CI of each mutation and write
    for eM in sP:
        mci = wilsonSI(sP[eM], tm)

        out.write(eM + "," + str(sP[eM]) + "," + str(mci[0]) + "," + str(mci[1]) + "\n")

    out.close()