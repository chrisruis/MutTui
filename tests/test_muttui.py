# test add_tree_node_labels
from MutTui.muttui import main
import sys
import os
import tempfile


def compare_files(file1, file2):
    with open(file1) as f1:
        with open(file2) as f2:
            assert f1.read() == f2.read()
    return


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
        # assert os.path.isfile(tmpdirname + "sequence_evolution_model.txt")

        # check output
        compare_files(tmpdirname + "DBS_label_A.csv", 
                      datafolder + 'expected_output/' + "DBS_label_A.csv")
        compare_files(tmpdirname + "all_included_double_substitutions.csv", 
                      datafolder + 'expected_output/' + "all_included_double_substitutions.csv")
        compare_files(tmpdirname + "all_included_mutations.csv", 
                      datafolder + 'expected_output/' + "all_included_mutations.csv")
        compare_files(tmpdirname + "ancestral_sequences.fasta", 
                      datafolder + 'expected_output/' + "ancestral_sequences.fasta")
        compare_files(tmpdirname + "annotated_tree.nexus", 
                      datafolder + 'expected_output/' + "annotated_tree.nexus")
        # compare_files(tmpdirname + "branch_mutations.txt", 
        #               datafolder + 'expected_output/' + "branch_mutations.txt")
        compare_files(tmpdirname + "gaps_to_N_alignment.fasta", 
                      datafolder + 'expected_output/' + "gaps_to_N_alignment.fasta")
        compare_files(tmpdirname + "mutation_types_label_A.csv", 
                      datafolder + 'expected_output/' + "mutation_types_label_A.csv")
        compare_files(tmpdirname + "mutational_spectrum_label_A.csv", 
                      datafolder + 'expected_output/' + "mutational_spectrum_label_A.csv")
        compare_files(tmpdirname + "mutations_not_included.csv", 
                      datafolder + 'expected_output/' + "mutations_not_included.csv")
        # compare_files(tmpdirname + "sequence_evolution_model.txt", 
        #              datafolder + 'expected_output/' + "sequence_evolution_model.txt")

    return

