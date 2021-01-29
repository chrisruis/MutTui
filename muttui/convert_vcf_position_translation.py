#Converts a VCF file created to a conversion between alignment position and genome position
#In output, column 1 is alignment position, column 2 is genome position

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-v",
                        "--vcf",
                        dest = "vcf",
                        required = True,
                        help = "VCF file to be converted")
    parser.add_argument("-f",
                        "--format",
                        dest = "format",
                        help = "The format of the VCF file to be converted. Options are:\n" + 
                        "gubbins\n" + 
                        "gubbins is the default",
                        default = "gubbins")
    parser.add_argument("-o",
                        "--out",
                        dest = "outfile",
                        required = True,
                        help = "Name of output file")
    args = parser.parse_args()

    #Import the vcf file
    vcf = open(args.vcf).readlines()

    outFile = open(args.outfile, "w")

    #Will be incremented with each alignment position
    iterator = 1

    #Check if the VCF file is in gubbins format
    if args.format == "gubbins":
        
        #Iterate through the VCF file and write the conversion between alignment and genome positions
        for line in vcf:
            if line.strip()[0] != "#":
                outFile.write(str(iterator) + "\t" + line.strip().split("\t")[1] + "\n")
                iterator += 1

    outFile.close()