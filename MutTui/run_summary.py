#Produces summaries of a MutTui run
#Run in the same directory as the MutTui output files
#Doesn't require any arguments as the all_included_mutations.csv and SBS spectra are identified automatically

from os.path import exists
from os import listdir

#Summarise all_included_mutations.csv
if exists("all_included_mutations.csv"):
    nm = 0
    mutations = dict()
    with open("all_included_mutations.csv") as f:
        next(f)
        for l in f:
            nm += 1

            if l.strip().split(",")[0] not in mutations:
                mutations[l.strip().split(",")[0]] = [l.strip().split(",")[1], 1]
            else:
                mutations[l.strip().split(",")[0]][1] += 1
    
    print("Total mutations in run:", nm)

    out = open("all_mutation_counts.csv", "w")
    out.write("Mutation_in_alignment,Mutation_in_genome,Number_of_mutations\n")

    for m in mutations:
        out.write(m + "," + mutations[m][0] + "," + str(mutations[m][1]) + "\n")

    out.close()
else:
    raise RuntimeError("all_included_mutations.csv file not present. This script should be run from the directory containing MutTui output files")

#Summarise SBS spectra
nf = 0
labelM = dict()
for f in listdir("."):
    if ("mutational_spectrum_" in f) and (".csv" in f) and ("_rescaled.csv" not in f):
        nf += 1
        label = f.replace("mutational_spectrum_label_", "").replace(".csv", "")
        labelM[label] = 0

        with open(f) as fn:
            next(fn)
            for l in fn:
                labelM[label] += int(l.strip().split(",")[1])

if nf == 0:
    raise RuntimeError("No SBS spectra found. This script should be run from the directory containing MutTui output files")
else:
    print("Number of analysed labels:", nf)
    print("Number of mutations per label:")
    for eL in labelM:
        print(eL + ":" + str(labelM[eL]))