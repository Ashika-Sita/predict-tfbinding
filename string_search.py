
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()
parser.add_argument('--f', '--input-fasta', dest='infile', type=str, help='input-fasta', metavar='input_fasta')
parser.add_argument('--m', '--motif', dest='mfile', type=str, help='motifs to search for', metavar='motif')
parser.add_argument('--o', '--outfile', dest='outfile', type=str, help='file name of output', metavar='outfile')
args = parser.parse_args()


def read_fasta(fastafile):

    """read in the input fasta file and return iterator over fasta records in file"""

    input_iterator = SeqIO.parse(fastafile, "fasta", alphabet = IUPAC.ambiguous_dna)

    return input_iterator

def search_fasta(input_iterator, motifile):

    """takes in fasta iterator and file of motifs;
    searches for motifs and
    returns dictionary of record:list of tuples of location, sequence for each record"""

    motif_dictionary = {}

    with open(motifile, 'r') as m:
        m_lines = m.readlines()

    instances = [Seq(str(line.strip('\n')), IUPAC.ambiguous_dna) for line in m_lines]
    motif_object = motifs.create(instances)

    for item in input_iterator:

        location_list = [(position, sequence) for position, sequence in motif_object.instances.search((item.seq).upper())]
        chrom_location = item.description.split(" ")[1]
        motif_dictionary[item.id+"--"+chrom_location] = location_list

    return motif_dictionary

def write_outfile(motif_dictionary, outfile):

    outfile = open(outfile, 'w')

    for item in motif_dictionary.keys():
        outfile.write(str(item.split("--")[0]) + "\t" + str(item.split("--")[1]) + "\t" + str(motif_dictionary[item]) + "\t" + str(len(motif_dictionary[item])) + "\n")
    outfile.close()

    return

fasta_files = read_fasta(args.infile)
searches = search_fasta(fasta_files, args.mfile)
write_outfile(searches, args.outfile)
