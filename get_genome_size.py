import sys
from Bio import SeqIO

usage = 'Usage: ' + sys.argv[0] + ' <fasta file> '

if len(sys.argv) < 2:
    print(usage)
    exit()

def main():
    fasta = sys.argv[1]
    data = SeqIO.parse(fasta, "fasta")
    for record in data:
        print(record.id)
        print(len(record.seq))

main()

