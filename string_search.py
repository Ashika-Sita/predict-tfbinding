from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC


def read_fasta(test_fasta):

    """read in the input fasta file and return iterator over fasta records in file"""

    input_iterator = SeqIO.parse(test_fasta, "fasta", alphabet = IUPAC.ambiguous_dna)

    return input_iterator

def search_fasta(string_file, input_iterator):

    """takes in fasta iterator and file of motifs;
    searches for motifs and
    returns dictionary of record:list of tuples of location, sequence for each record"""

    motif_dictionary = {}
    motif_list = []
    with open(string_file, 'r') as m:
        m_lines = m.readlines()

    instances = [Seq(str(line.strip('\n')), IUPAC.ambiguous_dna) for line in m_lines]
    motif_object = motifs.create(instances)
    motif_object.weblogo("./motif_search_outfiles/string_motif.png")

    for item in input_iterator:

        motif_list.append(item.id)

    return motif_list
