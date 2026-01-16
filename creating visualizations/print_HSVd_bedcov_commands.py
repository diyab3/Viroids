import sys

usage = 'Usage: ' + sys.argv[0] + ' <BED file> ' + '<BAM files>'

if len(sys.argv) < 3:
    print(usage)
    exit()

bed_file = sys.argv[1]
bam_files = sys.argv[2:] 
for bam in bam_files:
    command = f'samtools bedcov {bed_file} {bam}'
    print(command)


