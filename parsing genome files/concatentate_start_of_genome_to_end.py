import sys
from Bio import SeqIO

usage = 'Usage: ' + sys.argv[0] + ' <fastq file>'

if len(sys.argv) < 2:
    print(usage)
    exit()

fasta = sys.argv[1]

data = SeqIO.parse(fasta, "fasta")
for record in data:
    print(">" + record.id)
    new_end = record.seq[0:50]
    new_seq = record.seq[50:] + new_end

    print(new_seq)


