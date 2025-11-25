import sys
import os
usage = 'Usage: ' + sys.argv[0] + '<BAM files>'

if len(sys.argv) < 2:
    print(usage)
    exit()

bam_files = sys.argv[1:]
for bam in bam_files:
    basename,ext = os.path.splitext(bam)
    command = f'samtools depth {bam} > {basename}.depth'
    print(command)

