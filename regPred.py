"""
Input: 1) FASTA file of ChIP-seq identified binding sites
        2) File of string motif instances experimentally validated for TF binding
        3) FASTA file of sequences to search for TF binding


Output: All stored in folder motifSearchOutfiles
        1) String, PWM and TFFM logo PNGs
        2) Venn diagram of records predicted by all motif prediction algorithms

"""


from run_meme import *
from tffm_search import *
from string_search import *
from pwm_search import *
from draw_venn import *
from Bio import motifs
import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('-cf', '--chip-fasta', dest='chip_fasta', type=str, help='A FASTA file generated from ChIP-seq data for your transcription factor. Only uses the first 100 records, so order records by peak score')

parser.add_argument('-tf', '--test-fasta', dest='test_fasta', type=str, help='The FASTA file you wish to search for TFBS')
parser.add_argument('-sm','--string_motif', dest='string_motif', type=str, help='A text file with unambiguous DNA motifs; ideally these have been validated experimentally')


args = parser.parse_args()
subprocess.call("mkdir motif_search_outfiles", shell=True)
dictionary_for_venn = {}

run_meme(args.chip_fasta)

meme_file = "./temp/memeOut/meme.txt"

###run pwm search

pssm_list = pssm_search(meme_file, args.test_fasta)
dictionary_for_venn["pssm_search"] = pssm_list

###run string search

input_file = read_fasta(args.test_fasta)
string_list = search_fasta(args.string_motif, input_file)
dictionary_for_venn["string_search"] = string_list

###run tffm search

tffm_list = tffm_detailed(meme_file, args.chip_fasta, args.test_fasta)
dictionary_for_venn["tffm_search"] = tffm_list

###draw venn
draw_venn(dictionary_for_venn)

###delete temp

subprocess.call("rm -rf temp", shell=True)
