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

        location_list = [(position, sequence) for position, sequence in motif_object.instances.search((item.seq).upper())]
        chrom_location = item.description.split(" ")[1]
        motif_dictionary[item.id+"--"+chrom_location] = location_list
        motif_list.append(item.id)

    return motif_list

def write_outfile(motif_dictionary, outfile):

    outfile = open(outfile, 'w')

    for item in motif_dictionary.keys():
        outfile.write(str(item.split("--")[0]) + "\t" + str(item.split("--")[1]) + "\t" + str(motif_dictionary[item]) + "\t" + str(len(motif_dictionary[item])) + "\n")
    outfile.close()

    return
