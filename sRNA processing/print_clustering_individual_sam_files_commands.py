import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <sam files> '

if len(sys.argv) < 1:
    print(usage)
    exit()

sam_files = sys.argv[1:]
for sam in sam_files:
    sam_prefix, sam_suffix = os.path.splitext(sam)
    command = f'python print_read_clusters_v3.py {sam} > {sam_prefix}_clusters.sh'
    print(command)

