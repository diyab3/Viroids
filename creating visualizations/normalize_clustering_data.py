import sys
import pandas as pd
import os
usage = 'Usage: ' + sys.argv[0] + ' <clustering output files> '

if len(sys.argv) < 1:
    print(usage)
    exit()

clustering_files = sys.argv[1:]
total_reads = 1838360
HSVd_or_Cusa = "HSVd"

for cluster_file in clustering_files:
    with open(cluster_file, newline='') as tsvfile:
        df = pd.read_csv(cluster_file, sep="\t", header=None, names=['pos', 'strand', 'num reads', "5' het", 'seq'])
    
        df['RPM'] = (df['num reads'] / total_reads) * 1000000
        basename,ext = os.path.splitext(cluster_file)
        output_file = "normalized_" + basename + "_" + HSVd_or_Cusa + ".tsv"

        with open(output_file, 'w') as output:
            df.to_csv(output_file, sep="\t", index=False)
