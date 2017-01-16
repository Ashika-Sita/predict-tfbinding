from Bio import SeqIO
from Bio import Alphabet
from Bio import motifs

def pssm_search(meme_file, test_fasta):

    meme_handle = open(meme_file)
    meme_record = motifs.parse(meme_handle, "meme")
    meme_handle.close()
    meme_motif = meme_record[0]

    print meme_motif.degenerate_consensus

    meme_pssm = meme_motif.pssm
    meme_motif.weblogo("./motif_search_outfiles/pwm_motif.png")

    fasta_iterator = SeqIO.parse(test_fasta, "fasta", alphabet = meme_motif.alphabet)
    results = []

    for item in fasta_iterator:

        per_item_matches = []

        for pos, score in meme_pssm.search(item.seq, threshold = 3.0):
            per_item_matches.append(pos)

        if len(per_item_matches) > 0:
            results.append(item.id)

    return results
