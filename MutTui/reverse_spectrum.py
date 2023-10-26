#Reverses an RNA spectrum
#Can be used to generate a genomic strand spectrum from sequences stored as cRNA

import argparse
from reconstruct_spectrum import getRNADict, complement

if __name__ == "__main__":
    description = "Reverses an RNA mutational spectrum"
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument("-s",
                        "--spectrum",
                        dest = "spectrum",
                        required = True,
                        help = "Mutational spectrum to be reversed",
                        type = argparse.FileType("r"))
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Output csv to which the reversed spectrum will be written")
    
    args = parser.parse_args()

    #Generate empty spectrum to be filled
    s = getRNADict()

    #Iterate through the spectrum, reverse each mutation and add to s
    with open(args.spectrum.name) as f:
        next(f)
        for l in f:
            m = l.strip().split(",")[0]
            s[complement(m[6]) + complement(m[2]) + complement(m[4]) + complement(m[0])] = l.strip().split(",")[1]

    #Write the spectrum
    out = open(args.outFile, "w")
    out.write("Substitution,Number_of_mutations\n")
    for eM in s:
        out.write(eM[0] + "[" + eM[1] + ">" + eM[2] + "]" + eM[3] + "," + s[eM] + "\n")
    out.close()