import sys

usage = 'Usage: ' + sys.argv[0] + ' <lane info file> ' + '<bowtie2 index name>'

if len(sys.argv) != 3:
    print(usage)
    exit()

lane_info_file = sys.argv[1]
bowtie2_idx = sys.argv[2]

with open(lane_info_file) as F:
    for line in F:
        name, trim, suffix = line.split(".")
        command = f'hqsub "bowtie2 -x {bowtie2_idx} -U {name}.trim.fastq -S {name}.sam --un {name}_unaligned_reads.fastq" -q "*" -r {name}.log'
        
        print(command)

