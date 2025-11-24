import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <lane info file>'

if len(sys.argv) != 2:
    print(usage)
    exit()

lane_info_file = sys.argv[1]

with open(lane_info_file) as F:
    for line in F:
        sra_acc, name = line.strip().split()
        fastq_command = f'mv {sra_acc}.fastq {name}.fastq'
        print(fastq_command)
