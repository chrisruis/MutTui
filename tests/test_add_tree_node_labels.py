# test add_tree_node_labels
from muttui.add_tree_node_labels import main
import sys, os
import tempfile

def test_add_tree_node_labels(datafolder):

    with tempfile.TemporaryDirectory() as tmpdirname:

        print(os.listdir())
        # run add_tree_node_labels
        sys.argv = [
            "", "-t",  datafolder + "orygis_rooted.nwk", "-o",
            tmpdirname + "orygis_rooted_labelled.nwk"
        ]
        main()

        # check output
        with open(tmpdirname + "orygis_rooted_labelled.nwk", 'r') as infile:
            assert next(infile).strip()=='''((((((((((((SRR2100577.1:0.05954,\
SRR5642717:0.03823)Node12:0.00000,(ERR2659156:0.04048,ERR2659\
155:0.04448)Node13:0.00000)Node11:0.00141,SRR5642715:0.04889)\
Node10:0.00973,SRR2101222:0.05958)Node9:0.00303,(SRR5642712:0\
.05116,SRR5642711:0.04665)Node14:0.00147)Node8:0.00000,SRR210\
1291:0.04564)Node7:0.00359,SRR5642716:0.05632)Node6:0.00992,S\
RR2101161:0.06548)Node5:0.00078,ERR2659154:0.06354)Node4:0.00\
154,((SRR5642718:0.02170,SRR5642713:0.02653)Node16:0.06462,SR\
R5642714:0.07196)Node15:0.00000)Node3:0.00265,ERR2659153:0.09\
780)Node2:0.03955,SRR2101329:0.13735)Node1:0.00000;'''
    

    return