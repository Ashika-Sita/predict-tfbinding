from string_search import *
from pwm_search import *
from tffm_search import *
from draw_venn import *
dict_for_venn = {}
##string

fasta_files = read_fasta("test.fa")
string_searches = search_fasta(fasta_files, "meme_string_motif_test.txt")
dict_for_venn["string"] = string_searches

##pwm

pwm_searches = pssm_search("meme.txt", "test.fa")
dict_for_venn["pwm"] = pwm_searches

##tffm

tffm_searches = tffm_detailed("meme.txt", "train.fa", "test.fa")
dict_for_venn["tffm"] = tffm_searches


##create venn
draw_venn3(dict_for_venn)
