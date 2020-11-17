#Runs treetime ancestral reconstruction on an alignment and tree

import os
import subprocess

#Runs treetime ancestral reconstruction
def run_treetime(alignment, tree, output_dir, add_treetime_cmds):
    cmd = "treetime ancestral "

    #Check if a model has been specified, if not use --gtr infer
    if add_treetime_cmds == None:
        cmd += "--gtr infer"
    elif "--gtr" not in add_treetime_cmds:
        cmd += "--gtr infer"
    else:
        cmd += add_treetime_cmds
    
    cmd += " --aln " + alignment.name
    cmd += " --tree " + tree.name
    cmd += " --outdir " + output_dir

    subprocess.run(cmd, shell = True, check = True, executable = '/bin/sh')

    #Check if treetime completed successfully
    if (os.stat(output_dir + "ancestral_sequences.fasta").st_size == 0) or (os.stat(output_dir + "annotated_tree.nexus") == 0):
        raise Exception("treetime did not complete successfully")

#Runs treetime mugration
def run_treetime_mugration(tree, labels, output_dir):
    cmd = "treetime mugration --tree "

    cmd += tree
    cmd += " --states " + labels
    cmd += " --attribute group --confidence --outdir " + output_dir + "/mugration_out"

    subprocess.run(cmd, shell = True, check = True, executable = '/bin/sh')