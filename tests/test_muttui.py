# test add_tree_node_labels
from muttui.muttui import main
import sys
import os
import tempfile

def test_add_tree_node_labels(datafolder):

    with tempfile.TemporaryDirectory() as tmpdirname:

        # run add_tree_node_labels
        sys.argv = [
            "run", 
            "-a", datafolder + "orygis.filtered_polymorphic_sites.fasta", 
            "-t", datafolder + "orygis_rooted.nwk",
            "-r", datafolder + "SRR2100577.contigs_velvet.fa",
            "-c", datafolder + "conversion.txt",
            "-o", tmpdirname
        ]
        main()

        tmpdirname += '/'
        
        # check files are present
        assert os.path.isfile(tmpdirname + "DBS_label_A.csv")
        assert os.path.isfile(tmpdirname + "DBS_label_A.pdf")
        assert os.path.isfile(tmpdirname + "all_included_double_substitutions.csv")
        assert os.path.isfile(tmpdirname + "all_included_mutations.csv")
        assert os.path.isfile(tmpdirname + "ancestral_sequences.fasta")
        assert os.path.isfile(tmpdirname + "annotated_tree.nexus")
        assert os.path.isfile(tmpdirname + "branch_mutations.txt")
        assert os.path.isfile(tmpdirname + "gaps_to_N_alignment.fasta")
        assert os.path.isfile(tmpdirname + "mutation_types_label_A.csv")
        assert os.path.isfile(tmpdirname + "mutation_types_label_A.pdf")
        assert os.path.isfile(tmpdirname + "mutational_spectrum_label_A.csv")
        assert os.path.isfile(tmpdirname + "mutational_spectrum_label_A.pdf")
        assert os.path.isfile(tmpdirname + "mutations_not_included.csv")
        assert os.path.isfile(tmpdirname + "sequence_evolution_model.txt")

        # check output
        # TODO: add tests for output
    

    return

