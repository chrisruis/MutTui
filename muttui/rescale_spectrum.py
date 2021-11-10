#Rescales a mutational spectrum with a contexts file
#The relative height of each peak in the spectrum will be scaled by the proportion of the context in the given file

import argparse
from plot_spectrum import convertSpectrumDict
from reconstruct_spectrum import mutationTypeCount, getMutationDict, getRNADict
from calculate_contexts import *

#Removes the brackets and > from spectrum keys
def convertSpectrumKeys(spectrum):
    s = dict()

    for k in spectrum:
        s[k[0] + k[2] + k[4] + k[6]] = spectrum[k]

    return(s)

#Sums the contextual scalars for each mutaition type
def getMutationScalars(spectrum, contexts):
    #Mutation types as keys, total scalar as values
    scalars = dict()

    #Iterate through the mutations in the spectrum and add to scalars
    for m in spectrum:
        mutation = m[1:3]

        if mutation not in scalars:
            scalars[mutation] = float(0)
        
        scalars[mutation] += float(spectrum[m])/float(contexts[m[0] + m[3]])

    return(scalars)

#Rescales a mutation count by the number of contexts and, the total mutations of that type and the GC content
def rescaleMutation(x, c, s, tm, g):
    return(round(((x/c)/s) * tm * g))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--spectrum",
                        dest = "spectrum",
                        required = True,
                        help = "Mutational spectrum to be rescaled",
                        type = argparse.FileType("r"))
    
    sContext = parser.add_mutually_exclusive_group(required = True)
    sContext.add_argument("-c",
                        "--contexts",
                        dest = "contexts",
                        help = "File containing context proportions to use for the rescaling. Can be created using " + 
                        "calculate_contexts.py",
                        type = argparse.FileType("r"))
    sContext.add_argument("-r",
                        "--reference",
                        dest = "reference",
                        help = "Fasta format sequence to be examined. This should be a single sequence",
                        type = argparse.FileType("r"))
    
    parser.add_argument("-gc",
                        dest = "gc",
                        help = "GC content of the genome. Required if using -c, not required if using -r",
                        default = False)
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen, if specified the contexts will not include reverse complements",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        required = True,
                        help = "Output csv to which the rescaled spectrum will be written")
    
    args = parser.parse_args()

    #Check if no GC content is provided when providing contexts
    if args.contexts and not args.gc:
        raise RuntimeError("The proportion GC needs to be provided with -gc when using -c")

    #Import the spectrum and convert to a dictionary
    spectrum = convertSpectrumDict(args.spectrum)

    #Update spectrum keys
    spectrum = convertSpectrumKeys(spectrum)

    #Calculate the number of each mutation type
    mt = mutationTypeCount(spectrum, args.rna)

    if args.contexts:
        #Import the contexts and extract to dictionary
        contexts = dict()
        with open(args.contexts.name) as fileobject:
            next(fileobject)
            for line in fileobject:
                contexts[line.strip().split(",")[0]] = line.strip().split(",")[1]
        gc = float(args.gc)
        at = 1 - gc
    #Extract the contexts and GC content from the provided reference
    else:
        for record in SeqIO.parse(args.reference, "fasta"):
            contexts = calculateContexts(record.seq.upper(), args.rna)
            gcC, atC = getGC(record.seq.upper())
            gc = float(gcC)/(float(gcC) + float(atC))
            at = float(atC)/(float(gcC) + float(atC))
    
    #Sum the contextual scalars for each mutation
    mtScalars = getMutationScalars(spectrum, contexts)

    #Scalar for mutations starting with AT nucleotides
    atScalar = {"C":float(1), "G":float(1), "A":gc/at, "T":gc/at}

    #Will be filled with the rescaled spectrum
    if args.rna:
        rescaledSpectrum = getRNADict()
    else:
        rescaledSpectrum = getMutationDict()
    
    #Iterate through the contextual mutations, scale their counts and add to rescaledSpectrum
    for m in spectrum:
        rescaledSpectrum[m] = rescaleMutation(float(spectrum[m]), float(contexts[m[0] + m[3]]), mtScalars[m[1] + m[2]], float(mt[m[1] + m[2]]), atScalar[m[1]])
    
    outFile = open(args.outFile, "w")
    outFile.write("Substitution,Number_of_mutations\n")

    for eachMutation in rescaledSpectrum:
        outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(rescaledSpectrum[eachMutation]) + "\n")

    outFile.close()