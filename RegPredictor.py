import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from TFFM import tffm_module
from TFFM.constants import TFFM_KIND
from matplotlib import pyplot as plt
from matplotlib_venn import venn3_unweighted

class RegPred(object):


    def __init__(self, train_fasta, test_fasta, string_motif):

        self.train_fasta = train_fasta
        self.test_fasta = test_fasta
        self.string_motif = string_motif

    def run_meme(self):

        """Runs MEME on the train.fa file to identify top enriched
        motif for PWM search, and to initialize TFFM search"""

        counter =  0
        fasta_iterator = SeqIO.parse(self.train_fasta, "fasta")

        records = []

        for item in fasta_iterator:

            if counter < 100:

                records.append(item)
                counter += 1


        print len(records)
        SeqIO.write(records, "./temp-regPred/top100.fa", "fasta" )

        subprocess.call("meme ./temp-regPred/top100.fa -oc temp-regPred/memeOut -dna -mod anr -nmotifs 1 -minw 8 -maxw 20 -maxsize 1000000000, -revcomp", shell=True)



        return

    def create_fa_iter(self):

        """Creates an iterator object over records in the test.fa file
        for use by motif search modules"""

        input_iterator = SeqIO.parse(self.test_fasta, "fasta", alphabet = IUPAC.ambiguous_dna)

        print "FASTA iterator created"

        return input_iterator

    def string_search(self, input_iterator):

        """Searches for exact matches to user inputted experimentally
        validated motifs"""

        print "Starting String Search"

        string_list = []
        outfile = open("./motif_search_outfiles/string_ids.txt", 'w')

        with open(self.string_motif, 'r') as m:
            m_lines = m.readlines()

        instances = [Seq(str(line.strip('\n')), IUPAC.ambiguous_dna) for line in m_lines]
        motif_object = motifs.create(instances)
        print "String Motif:", motif_object.degenerate_consensus
        motif_object.weblogo("./motif_search_outfiles/string_motif.png")

        print "Parsing for string matches"

        for item in input_iterator:
            location_list = [(position, sequence) for position, sequence in motif_object.instances.search((item.seq).upper())]
            if len(location_list) > 0:
                string_list.append(item.id)

        for item in string_list:
            outfile.write(str(item))
            outfile.write("\n")
        outfile.close()

        print "Completed string search"

        return string_list

    def pssm_search(self, input_iterator):

        """Searches for PSSM matches above threshold 3.0 in test sequences"""
        """Creates an iterator object over records in the test.fa file
        for use by motif search modules"""

        #input_iterator = SeqIO.parse(self.test_fasta, "fasta", alphabet = IUPAC.ambiguous_dna)

        #print "FASTA iterator created"

        print "Starting PSSM Search"

        pssm_list = []
        outfile = open("./motif_search_outfiles/pssm_ids.txt", 'w')


        meme_handle = open("./temp-regPred/memeOut/meme.txt", 'r')
        meme_record = motifs.parse(meme_handle, "meme")
        meme_handle.close()
        meme_motif = meme_record[0]
        meme_pssm = meme_motif.pssm
        meme_motif.weblogo("./motif_search_outfiles/pwm_motif.png")

        print "PSSM Motif:", meme_motif.degenerate_consensus

        print "Parsing for PSSM matches"

        for item in input_iterator:

            item.seq.alphabet = meme_motif.alphabet
            per_item_matches = []

            for pos, score in meme_pssm.search(item.seq, threshold = 3.0):

                per_item_matches.append(pos)

            if len(per_item_matches) > 0:
                pssm_list.append(item.id)

        #print pssm_list

        for item in pssm_list:
            outfile.write(str(item))
            outfile.write("\n")

        outfile.close()

        print "Completed PSSM search"

        return pssm_list

    def tffm_search(self):

        print "Starting TFFM search"

        tffm_list = []
        outfile = open("./motif_search_outfiles/tffm_ids.txt", 'w')

        tffm_detailed = tffm_module.tffm_from_meme("./temp-regPred/memeOut/meme.txt", TFFM_KIND.DETAILED)
        tffm_detailed.write("./temp-regPred/tffm_detailed_initial.xml")
        tffm_detailed.train(self.train_fasta)
        tffm_detailed.write("./temp-regPred/tffm_detailed.xml")
        out = open("./motif_search_outfiles/tffm_detailed_summary_logo.svg", "w")
        tffm_detailed.print_summary_logo(out)
        out.close()
        out = open("./motif_search_outfiles/tffm_detailed_dense_logo.svg", "w")
        tffm_detailed.print_dense_logo(out)
        out.close()

        tffm_detailed = tffm_module.tffm_from_xml("./temp-regPred/tffm_detailed.xml",
            TFFM_KIND.DETAILED)

        print "Parsing for TFFM matches"

        for hit in tffm_detailed.scan_sequences(self.test_fasta, only_best=True):
            tffm_list.append(str(hit).split('\t')[0])

        for item in tffm_list:
            outfile.write(str(item))
            outfile.write("\n")

        outfile.close()

        print "Completed TFFM search"

        return tffm_list

    def find_common(self, string_list, pssm_list, tffm_list):

        outfile = open("./motif_search_outfiles/common_to_all.txt", 'w')
        common_to_three = set.intersection(set(string_list), set(pssm_list), set(tffm_list))

        for item in common_to_three:
            outfile.write(str(item))
            outfile.write("\n")


        outfile.close()

        return

    def draw_venn(self, string_list, pssm_list, tffm_list):

        venn3_unweighted([set(string_list), set(pssm_list), set(tffm_list)], ["String-Search", "PSSM-Search", "TFFM-Search"])
        plt.savefig("./motif_search_outfiles/venn_motifSearches.png")

        return

subprocess.call("rm -rf temp-regPred", shell=True)
