#Calculates the number of each possible context from a genome sequence

import argparse
from Bio import SeqIO
from collections import OrderedDict
from plot_spectrum import convertSpectrumProportions

#Creates a dictionary of contexts
def getContextDict():
    contexts = OrderedDict()
    contexts["AA"] = 0
    contexts["AC"] = 0
    contexts["AG"] = 0
    contexts["AT"] = 0
    contexts["CA"] = 0
    contexts["CC"] = 0
    contexts["CG"] = 0
    contexts["CT"] = 0
    contexts["GA"] = 0
    contexts["GC"] = 0
    contexts["GG"] = 0
    contexts["GT"] = 0
    contexts["TA"] = 0
    contexts["TC"] = 0
    contexts["TG"] = 0
    contexts["TT"] = 0

    return(contexts)

#Calculates the proportion of each context in a given sequence
def calculateContexts(sequence, rna):
    contexts = getContextDict()

    #Iterate through the sequence and add the contexts to contexts
    for s in range(len(sequence) - 2):
        contexts[sequence[s] + sequence[s+2]] += 1
    
    #Add the reverse contexts if using DNA
    if not rna:
        withRC = getContextDict()
        withRC["AA"] = contexts["AA"] + contexts["TT"]
        withRC["AC"] = contexts["AC"] + contexts["GT"]
        withRC["AG"] = contexts["AG"] + contexts["CT"]
        withRC["AT"] = contexts["AT"] + contexts["AT"]
        withRC["CA"] = contexts["CA"] + contexts["TG"]
        withRC["CC"] = contexts["CC"] + contexts["GG"]
        withRC["CG"] = contexts["CG"] + contexts["CG"]
        withRC["CT"] = contexts["CT"] + contexts["AG"]
        withRC["GA"] = contexts["GA"] + contexts["TC"]
        withRC["GC"] = contexts["GC"] + contexts["GC"]
        withRC["GG"] = contexts["GG"] + contexts["CC"]
        withRC["GT"] = contexts["GT"] + contexts["AC"]
        withRC["TA"] = contexts["TA"] + contexts["TA"]
        withRC["TC"] = contexts["TC"] + contexts["GA"]
        withRC["TG"] = contexts["TG"] + contexts["CA"]
        withRC["TT"] = contexts["TT"] + contexts["AA"]
        sP = convertSpectrumProportions(withRC)
    else:
        sP = convertSpectrumProportions(contexts)
    
    return(sP)
    
    #print(contexts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s",
                        "--sequence",
                        dest = "sequence",
                        required = True,
                        help = "Fasta format sequence to be examined. This should be a single sequence",
                        type = argparse.FileType("r"))
    parser.add_argument("--rna",
                        dest = "rna",
                        help = "Specify if using an RNA pathogen, if specified the contexts will not include reverse complements",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "--outfile",
                        dest = "outFile",
                        help = "Output csv to which contexts will be written")
    
    args = parser.parse_args()

    #Extract all of the contexts from the given sequence
    for record in SeqIO.parse(args.sequence, "fasta"):
        contexts = calculateContexts(record.seq.upper(), args.rna)
    
    outFile = open(args.outFile, "w")
    outFile.write("Context,Proportion\n")

    #Write the contexts
    for c in contexts:
        outFile.write(c + "," + str(contexts[c]) + "\n")

    outFile.close()