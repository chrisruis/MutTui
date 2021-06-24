#Calculates a mutational spectrum from variant input data
#Can take VCF input or a file containing variants
#This variant file can be in any format as long as it has columns:
#position - containing genome position of mutation
#reference - reference base at position
#variant - variant base at position
#A VCF file is the default input, to use a variant file use --variant

import argparse
import pandas as pd
from Bio import SeqIO
import sys
from reconstruct_spectrum import getMutationDict, getContext, complement
from plot_spectrum import convertSpectrumFormat, plotSpectrumFromDict

#Extracts variants from a variant file, returns a list of lists with
#each list containing position, reference, variant
def extractFromVariantFile(vf):
    #Import variants
    variants = pd.read_csv(vf.name)

    #Check if the required columns are present
    if ("position" not in variants.columns) or ("reference" not in variants.columns) or ("variant" not in variants.columns):
        print("Variant file needs to contain columns position, reference and variant. Exiting")
        sys.exit()

    #Empty list to be filled
    vL = list()

    #Iterate through the variants and add to vL
    for i in range(variants.shape[0]):
        vL.append([variants["reference"][i], variants["variant"][i], variants["position"][i]])
    
    return(vL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-v",
                        "--variant_file",
                        dest = "v_file",
                        required = True,
                        help = "Input file containing variants. Expects a VCF file by default " +
                        "or a variant file if --variant is used",
                        type = argparse.FileType("r"))
    parser.add_argument("-r",
                        "--reference",
                        dest = "reference",
                        required = True,
                        help = "Fasta reference file. This is not updated so mutational contexts will be " + 
                        "extracted from this as it is provided",
                        type = argparse.FileType("r"))
    parser.add_argument("--variant",
                        dest = "variant",
                        help = "Use this to supply a variant file rather than a VCF. This variant file " + 
                        "can be in any format as long as it has columns: position (containing the genome position " + 
                        "of the variant), reference (containing the reference nucleotide at the position) and variant " + 
                        "(containing the variant nucleotide at the position)",
                        action = "store_true",
                        default = False)
    parser.add_argument("-o",
                        "-output_prefix",
                        dest = "out",
                        required = True,
                        help = "Prefix for output files. By default, these will be saved in the current directory")
    
    args = parser.parse_args()

    #Extract variants from input file
    if args.variant:
        variants = extractFromVariantFile(args.v_file)
    
    #Import reference and extract to string
    reference = str(SeqIO.read(args.reference.name, "fasta").seq).upper()

    outMutationsNotUsed = open(args.out + "_mutations_not_included.csv", "w")
    outMutationsNotUsed.write("Mutation_in_genome,Reason_not_included\n")

    outAllMutations = open(args.out + "_all_included_mutations.csv", "w")
    outAllMutations.write("Mutation_in_genome,Substitution\n")
    
    #Empty spectrum
    spectrum = getMutationDict()

    #Nucleotides, used to check if mutation involves 2 nucleotides
    nucleotides = ["A","C","G","T"]

    #Iterate through the variants, get their genomic context, check if they involve all nucleotides and add to spectrum
    for mutation in variants:
        if (mutation[0] not in nucleotides) or (mutation[1] not in nucleotides):
            outMutationsNotUsed.write(mutation[0] + str(mutation[2]) + mutation[1] + ",mutation_does_not_involve_2_nucleotides\n")
        else:
            mutationContext = getContext(mutation, reference)

            #Check if the mutation, upstream or downstream nucleotides are not A, C, G or T
            if (mutationContext[0] not in nucleotides) or (mutationContext[1] not in nucleotides):
                outMutationsNotUsed.write(mutation[0] + str(mutation[2]) + mutation[1] + ",surrounding_position_not_nucleotide\n")
            else:
                #This will be true for all RNA mutations and half of DNA mutations
                if (mutationContext[0] + mutation[0] + mutation[1] + mutationContext[1]) in spectrum:
                    spectrum[mutationContext[0] + mutation[0] + mutation[1] + mutationContext[1]] += 1
                    outAllMutations.write(mutation[0] + str(mutation[2]) + mutation[1] + "," + mutationContext[0] + "[" + mutation[0] + ">" + mutation[1] + "]" + mutationContext[1] + "\n")
                #Add to the corresponding complement
                else:
                    spectrum[complement(mutationContext[1]) + complement(mutation[0]) + complement(mutation[1]) + complement(mutationContext[0])] += 1
                    outAllMutations.write(complement(mutation[0]) + str(mutation[2]) + complement(mutation[1]) + "," + complement(mutationContext[1]) + "[" + complement(mutation[0]) + ">" + complement(mutation[1]) + "]" + complement(mutationContext[0]) + "\n")
    
    #Write the spectrum
    outFile = open(args.out + "_mutational_spectrum.csv", "w")
    outFile.write("Substitution,Number_of_mutations\n")
    for eachMutation in spectrum:
        outFile.write(eachMutation[0] + "[" + eachMutation[1] + ">" + eachMutation[2] + "]" + eachMutation[3] + "," + str(spectrum[eachMutation]) + "\n")
    outFile.close()

    #Plot the spectrum
    outSpectrum = open(args.out + "_mutational_spectrum.pdf", "w")
    spectrumFormat = convertSpectrumFormat(spectrum)
    plotSpectrumFromDict(spectrumFormat, outSpectrum)
    outSpectrum.close()