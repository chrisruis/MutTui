#Rescales a mutational spectrum with a contexts file
#The relative height of each peak in the spectrum will be scaled by the proportion of the context in the given file

import argparse
from collections import OrderedDict, Counter
from Bio import SeqIO
from plot_spectrum import convertSpectrumDict
from reconstruct_spectrum import getMutationDict, getRNADict, getDoubleSubstitutionDict, complement

#Removes the brackets and > from spectrum keys
def convertSpectrumKeys(spectrum):
    s = dict()

    for k in spectrum:
        s[k[0] + k[2] + k[4] + k[6]] = spectrum[k]

    return(s)

#Removes the > from keys in a double substitution spectrum
def convertDoubleKeys(spectrum):
    s = dict()

    for k in spectrum:
        s[k[0] + k[1] + k[3] + k[4]] = spectrum[k]
    
    return(s)

#Creates an empty dictionary of triplets
def getTripletDict():
    tripletDict = OrderedDict()
    tripletDict["ACA"] = 0
    tripletDict["ACC"] = 0
    tripletDict["ACG"] = 0
    tripletDict["ACT"] = 0
    tripletDict["CCA"] = 0
    tripletDict["CCC"] = 0
    tripletDict["CCG"] = 0
    tripletDict["CCT"] = 0
    tripletDict["GCA"] = 0
    tripletDict["GCC"] = 0
    tripletDict["GCG"] = 0
    tripletDict["GCT"] = 0
    tripletDict["TCA"] = 0
    tripletDict["TCC"] = 0
    tripletDict["TCG"] = 0
    tripletDict["TCT"] = 0
    tripletDict["ATA"] = 0
    tripletDict["ATC"] = 0
    tripletDict["ATG"] = 0
    tripletDict["ATT"] = 0
    tripletDict["CTA"] = 0
    tripletDict["CTC"] = 0
    tripletDict["CTG"] = 0
    tripletDict["CTT"] = 0
    tripletDict["GTA"] = 0
    tripletDict["GTC"] = 0
    tripletDict["GTG"] = 0
    tripletDict["GTT"] = 0
    tripletDict["TTA"] = 0
    tripletDict["TTC"] = 0
    tripletDict["TTG"] = 0
    tripletDict["TTT"] = 0

    return(tripletDict)

#Creates an empty dictionary of dinucleotides
def getDinucleotideDict():
    dnDict = OrderedDict()
    dnDict["AA"] = 0
    dnDict["AC"] = 0
    dnDict["AG"] = 0
    dnDict["AT"] = 0
    dnDict["CA"] = 0
    dnDict["CC"] = 0
    dnDict["CG"] = 0
    dnDict["GA"] = 0
    dnDict["GC"] = 0
    dnDict["TA"] = 0

    return(dnDict)

#Calculates the number of each triplet in a given sequence
def calculateContexts(sequence, rna):
    tripletDict = getTripletDict()

    nucleotides = ["A", "C", "G", "T"]

    for i in range(len(sequence) - 2):
        t = sequence[i:(i+3)]
        if (t[0] in nucleotides) and (t[1] in nucleotides) and (t[2] in nucleotides):
            if t in tripletDict:
                tripletDict[t] += 1
            else:
                tripletDict[complement(t[2]) + complement(t[1]) + complement(t[0])] += 1
    
    #Add the reverse complement contexts from the reverse strand if the sequence is DNA
    if not rna:
        for c in tripletDict:
            tripletDict[c] = tripletDict[c] * 2
    
    return(tripletDict)

#Calculates the number of each dinucleotide in a given sequence
def calculateDinucleotides(sequence, rna):
    dnDict = getDinucleotideDict()

    for i in range(len(sequence) - 1):
        if (sequence[i] + sequence[i + 1]) in dnDict:
            dnDict[sequence[i] + sequence[i + 1]] += 1
        else:
            dnDict[complement(sequence[i + 1]) + complement(sequence[i])] += 1
    
    #Add the reverse complement contexts from the reverse strand if the sequence is DNA
    if not rna:
        for c in dnDict:
            dnDict[c] = dnDict[c] * 2
    
    return(dnDict)

#Rescales a SBS mutational spectrum
def rescaleSBS(spectrum, reference, scalar, rna):
    #Update spectrum keys
    spectrum = convertSpectrumKeys(spectrum)

    #Extract all triplet contexts from the reference
    for record in SeqIO.parse(reference, "fasta"):
        contexts = calculateContexts(record.seq.upper(), rna)
    
    #Will be filled with the rescaled spectrum
    if args.rna:
        rescaledSpectrum = getRNADict()
    else:
        rescaledSpectrum = getMutationDict()
    
    #Iterate through the contextual mutations, scale their counts and add to rescaledSpectrum
    for m in spectrum:
        rescaledSpectrum[m] = round((spectrum[m]/contexts[m[0] + m[1] + m[3]]) * scalar)

    return(rescaledSpectrum)

#Rescales a mutation type spectrum
def rescaleMT(spectrum, reference, scalar, rna):
    #Count the number of each nucleotide in the refernece, combine AT and GC if not rna
    for record in SeqIO.parse(reference, "fasta"):
        c = Counter(record.seq.upper())
    
    conversion = dict()
    if rna:
        conversion["A"] = c["A"]
        conversion["C"] = c["C"]
        conversion["G"] = c["G"]
        conversion["T"] = c["T"]
    else:
        conversion["A"] = c["A"] + c["T"]
        conversion["C"] = c["C"] + c["G"]
        conversion["G"] = c["C"] + c["G"]
        conversion["T"] = c["A"] + c["T"]
    
    #Empty mutation type dicts
    if rna:
        rescaledSpectrum = {"AC": 0, "AG": 0, "AT": 0, "CA": 0, "CG": 0, "CT": 0, "GA": 0, "GC": 0, "GT": 0, "TA": 0, "TC": 0, "TG": 0}
    else:
        rescaledSpectrum = {"CA": 0, "CG": 0, "CT": 0, "TA": 0, "TC": 0, "TG": 0}
    
    #Rescale the mutations based on their original nucleotide
    for m in spectrum:
        rescaledSpectrum[m] = round((spectrum[m]/conversion[m[0]]) * scalar)
    
    return(rescaledSpectrum)

#Rescales a double substitution spectrum
def rescaleDouble(spectrum, reference, scalar, rna):
    #Convert spectrum keys
    spectrum = convertDoubleKeys(spectrum)

    #Count dinucleotide pairs in reference
    for sequence in SeqIO.parse(reference, "fasta"):
        dn = calculateDinucleotides(sequence.seq, rna)

    #Empty DBS dict
    rescaledSpectrum = getDoubleSubstitutionDict()

    #Iterate through the mutations and rescale by their starting context availability
    for m in spectrum:
        rescaledSpectrum[m] = round((spectrum[m]/dn[m[0] + m[1]]) * scalar)
    
    return(rescaledSpectrum)

if __name__ == "__main__":
    description = "Rescales a mutational spectrum based on triplet availability in a given reference"
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument("-s",
                        "--spectrum",
                        dest = "spectrum",
                        required = True,
                        help = "Mutational spectrum to be rescaled",
                        type = argparse.FileType("r"))
    parser.add_argument("-r",
                        "--reference",
                        dest = "reference",
                        help = "Fasta format sequence to be examined. This should be a single sequence",
                        type = argparse.FileType("r"))
    parser.add_argument("-f",
                        "--scale_factor",
                        dest = "scale_factor",
                        help = "The scalar each mutation/context will be multiplied by to obtain a rescaled mutation count, default 1,000,000. " + 
                        "This may need to be increased if the number of mutations is very small compared to the genome length",
                        default = "1000000")
    parser.add_argument("--mutation_type",
                        dest = "mt",
                        help = "Specify if rescaling a spectrum of the 6 mutation types rather than a 96 bar spectrum of contextual mutations. The " + 
                        "input spectrum will be rescaled by the number of each starting nucleotide in the genome",
                        action = "store_true",
                        default = False)
    parser.add_argument("--double",
                        dest = "double",
                        help = "Specify if rescaling a double substitution spectrum rather than a 96 bar spectrum of contextual mutations. The " + 
                        "input spectrum will be rescaled by the number of each dinucleotide pair in the genome",
                        action = "store_true",
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

    #Import the spectrum and convert to a dictionary
    spectrum = convertSpectrumDict(args.spectrum)

    #The scalar each mutation/context will be multiplied by to obtain a rescaled mutation count
    scalar = float(args.scale_factor)

    outFile = open(args.outFile, "w")
    outFile.write("Substitution,Number_of_mutations\n")

    if args.mt:
        #Rescale the mutation type spectrum
        rescaledSpectrum = rescaleMT(spectrum, args.reference, scalar, args.rna)
        for eachMutation in rescaledSpectrum:
            outFile.write(eachMutation + "," + str(rescaledSpectrum[eachMutation]) + "\n")
    elif args.double:
        rescaledSpectrum = rescaleDouble(spectrum, args.reference, scalar, args.rna)
        for eachMutation in rescaledSpectrum:
            outFile.write(eachMutation[:2] + ">" + eachMutation[2:] + "," + str(rescaledSpectrum[eachMutation]) + "\n")
    else:
        #Rescale the SBS spectrum
        rescaledSpectrum = rescaleSBS(spectrum, args.reference, scalar, args.rna)
        #Write the rescaled spectrum
        for eachMutation in rescaledSpectrum:
            outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(rescaledSpectrum[eachMutation]) + "\n")

    outFile.close()