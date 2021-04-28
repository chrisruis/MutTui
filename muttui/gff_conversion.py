#Functions to work with GFF files

from io import StringIO
from Bio import SeqIO
import gffutils as gff

#Clean other "##" starting lines from gff file, as it confuses parsers
#Taken from prokka.py in Panaroo
def clean_gff_string(gff_string):
    splitlines = gff_string.splitlines()
    lines_to_delete = []
    for index in range(len(splitlines)):
        if '##sequence-region' in splitlines[index]:
            lines_to_delete.append(index)
    for index in sorted(lines_to_delete, reverse=True):
        del splitlines[index]
    cleaned_gff = "\n".join(splitlines)
    return cleaned_gff

#Takes a GFF file and returns a list of lists
#Each entry in list is a gene with 4 components - gene name, gene start, gene end, strand
def convertGFF(gff_file_name):
    gff_file = open(gff_file_name, "r")

    #Open file, split into genes and sequence
    lines = gff_file.read().replace(",", "")
    split = lines.split("##FASTA")

    if len(split) != 2:
        print("Problem reading GFF file", gff_file.name)
        raise RuntimeError("Error reading GFF file")
    
    #Convert GFF sequence to fasta format
    with StringIO(split[1]) as temp_fasta:
        sequences = list(SeqIO.parse(temp_fasta, "fasta"))
    
    parsed_gff = gff.create_db(clean_gff_string(split[0]),
                                dbfn = ":memory:",
                                force = True,
                                keep_order = True,
                                from_string = True)
    
    geneCoordinates = dict()
    geneDict = dict()
    
    for entry in parsed_gff.all_features(featuretype = ()):
        if "CDS" not in entry.featuretype:
            continue
        else:
            geneCoordinates[entry.id] = [entry.start, entry.stop, entry.strand]
            for position in range(entry.start, (entry.stop + 1)):
                #If the position is in overlapping genes, all genes will be included separated by ____
                if position in geneDict:
                    geneDict[position] = geneDict[position] + "____" + entry.id
                else:
                    geneDict[position] = entry.id

    return(geneCoordinates, geneDict)