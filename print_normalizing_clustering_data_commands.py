import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <clustering output files> '

if len(sys.argv) < 1:
    print(usage)
    exit()

cluster_files = sys.argv[1:]
aligned_to_Cusa_or_HSVd = HSVd

for cluster_file in cluster_files:
    base_name, ext = os.path.splitext(cluster_file)
    command = f'python normalize_clustering_data.py {cluster_file} > normalized_{base_name}_{aligned_to_Cusa_or_HSVd}{ext}'
    print(command)

