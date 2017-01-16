"""Identifies fasta records with specified TFBS using the TFFM module from the Wasserman lab """

from TFFM import tffm_module
import sys
from TFFM import *
from TFFM.constants import TFFM_KIND

def tffm_first(motif_file, training_fasta, test_fasta):

    tffm_first_order = tffm_module.tffm_from_meme(motif_file, TFFM_KIND)
    tffm_first_order.write("tffm_first_order_initial.xml")
    tffm_first_order.train(training_fasta)
    tffm_first_order.write("tffm_first_order.xml")
    out = open("tffm_first_order_summary_logo.svg", 'w')
    tffm_first_order.print_summary_logo(out)
    out.close()
    out = open("tffm_first_order_dense_logo", 'w')
    tffm_first_order.print_dense_logo(out)
    out.close()

    tffm_first_order = tffm_module.tffm_from_xml("tffm_first_order.xml",
        TFFM_KIND.FIRST_ORDER)
    print "1st-order all"
    for hit in tffm_first_order.scan_sequences(test_fasta):
        if hit:
            print hit

    print "1st-order best"
    for hit in tffm_first_order.scan_sequences(test_fasta, only_best=True):
        print hit

    return

def tffm_detailed(meme_file, training_fasta, test_fasta):
    tffm_detailed = tffm_module.tffm_from_meme(meme_file, TFFM_KIND.DETAILED)
    tffm_detailed.write("./temp/tffm_detailed_initial.xml")
    tffm_detailed.train(training_fasta)
    tffm_detailed.write("./temp/tffm_detailed.xml")
    out = open("./motif_search_outfiles/tffm_detailed_summary_logo.svg", "w")
    tffm_detailed.print_summary_logo(out)
    out.close()
    out = open("./motif_search_outfiles/tffm_detailed_dense_logo.svg", "w")
    tffm_detailed.print_dense_logo(out)
    out.close()
    hits = []
    tffm_detailed = tffm_module.tffm_from_xml("./temp/tffm_detailed.xml",
        TFFM_KIND.DETAILED)

    for hit in tffm_detailed.scan_sequences(test_fasta, only_best=True):

        hits.append(str(hit).split('\t')[0])

    return hits
