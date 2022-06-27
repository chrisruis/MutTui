#Creates a spectrum containing only the mutation types that are elevated in one spectrum relative to another
#Used to create a condition-specific spectrum to examine the signatures that are more active in one condition relative to another
#The generated spectrum will contain the number of mutations provided with -n, default 10000

import argparse
from plot_spectrum import convertSpectrumDict
from compare_spectra import convertSpectrumProportions

#Calculates the difference in mutation proportions between 2 spectra
def sDifference(s1, s2, mutations):
    sDiff = dict()

    for m in mutations:
        sDiff[m] = s1[m] - s2[m]
    
    return(sDiff)

#Calculates the proportion of mutations to be assigned to each mutation type in a comparison spectrum
#Takes a comparison spectrum from sDifference
#Keeps mutations that are >0 in this comparison spectrum and converts each mutation to proportion of that
#mutation in the enriched mutations
def getMutationProportions(s, mutations):
    #Sum of the proportions that are >0
    totMut = float(0)

    #Iterate through the mutations, check if they are >0 in s, add to totMut if they are, change to 0 if they are not
    for m in mutations:
        if s[m] > 0.0:
            totMut += s[m]
        else:
            s[m] = 0.0
    
    #Iterate through the mutations and change to their scaled value
    #Will be 0 if they were <=0 in the original spectrum
    #Will be scaled to their original value/total proportion of mutations >0 if they were originally >0
    for mut in mutations:
        s[mut] = s[mut]/totMut
    
    return(s)

if __name__ == "__main__":
    description = "Creates a spectrum containing mutations elevated in one spectrum over another, in proportion to their elevation"
    parser = argparse.ArgumentParser(description = description)

    parser.add_argument("-s1",
                        "--spectrum1",
                        dest = "spectrum1",
                        required = True,
                        help = "Spectrum csv from MutTui. Mutation types enriched in this spectrum above " + 
                        "the spectrum provided with -s2 will be included in the output spectrum",
                        type = argparse.FileType("r"))
    parser.add_argument("-s2",
                        "--spectrum2",
                        dest = "spectrum2",
                        required = True,
                        help = "Spectrum csv from MutTui. This spectrum is used as the background. " + 
                        "Mutations in the spectrum provided with -s1 that are above this spectrum will be " + 
                        "included in the output spectrum",
                        type = argparse.FileType("r"))
    parser.add_argument("-n",
                        "--number_mutations",
                        dest = "n_mut",
                        help = "The number of mutations in the output spectrum, default is 10,000",
                        default = "10000")
    parser.add_argument("-o",
                        "--output_spectrum",
                        dest = "out_file",
                        required = True,
                        help = "File to which the mutations enriched in spectrum -s1 relative to spectrum " + 
                        "-s2 will be written")
    
    args = parser.parse_args()

    #Import the spectra and convert them to dictionaries
    s1 = convertSpectrumDict(args.spectrum1)
    s2 = convertSpectrumDict(args.spectrum2)
    #Convert to dictionaries of proportions
    sp1 = convertSpectrumProportions(s1)
    sp2 = convertSpectrumProportions(s2)

    #The number of mutations in the enriched spectrum
    nMut = float(args.n_mut)

    outFile = open(args.out_file, "w")
    outFile.write("Substitution,Sample1\n")

    #The 96 DNA mutations as a list, used so they are always in the same order
    mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    #Subtract spectrum 2 from spectrum 1
    sDiff = sDifference(sp1, sp2, mutations)

    #Convert the mutation proportions so mutations more common in -s2 are set to 0 and
    #mutations more common in -s1 are set to their relative proportion in the mutations that are enriched in -s1
    sEnriched = getMutationProportions(sDiff, mutations)

    #Write the mutations
    for mutation in mutations:
        outFile.write(mutation + "," + str(int(round((sEnriched[mutation] * nMut),0))) + "\n")

    outFile.close()