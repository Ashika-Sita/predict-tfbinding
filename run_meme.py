import subprocess
from Bio import SeqIO
from Bio import SeqRecord


subprocess.call("mkdir temp", shell=True)

def run_meme(fasta_chip):
    counter = 0

    while counter < 100:
        fasta_iterator = SeqIO.parse(fasta_chip, "fasta")
        counter +=1

    for item in fasta_iterator:
        SeqIO.write(item, "./temp/fasta_chip_top100.fa", "fasta")

    subprocess.call("meme ./temp/fasta_chip_top100.fa -oc ./temp/memeOut -dna -mod anr -nmotifs 1 -minw 6 -maxw 10 -maxsize 10000000", shell=True)

    return
