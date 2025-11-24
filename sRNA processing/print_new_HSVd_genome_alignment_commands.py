import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <bowtie2 index name> ' + '<fastq files>'

if len(sys.argv) < 3:
    print(usage)
    exit()

fasta_files = sys.argv[2:]
bowtie2_idx = sys.argv[1]

for fasta in fasta_files:
    basename = os.path.splitext(fasta)[0]
       
    command = f'hqsub "bowtie2 -x {bowtie2_idx} -U {fasta} -S {basename}.sam" -q "*" -r {basename}.log'
        
    print(command)
